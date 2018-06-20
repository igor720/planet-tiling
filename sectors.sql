-- Description:  SQL instructions for creating ocean and
--               land sectors (tiles, tesselation)
--               on spherical planet survice
-- license: MIT
-- license-file: LICENSE
--
-- verified on the configuration:
--   PostgreSQL 10.4
--   POSTGIS="2.4.4 r16526" PGSQL="100" GEOS="3.5.1-CAPI-1.9.1 r4246"
--     PROJ="Rel. 4.9.2, 08 September 2015" GDAL="GDAL 1.11.3,
--     released 2015/09/16" LIBXML="2.9.3" LIBJSON="0.11.99"
--     LIBPROTOBUF="1.2.1" RASTER
--
-- Instructions:
--
-- The following commands suggest default PostgreSQL+PostGIS
-- installation. 'riversz' and 'lands' shapefiles should be fetched
-- from forgedmaps.com or they must contain data with similar structure
--
-- Prepare import sql files (we use 'map' schema for planet map data):
-- shp2pgsql -s 3785 -d -D -I -t 2D -N skip riversz.shp map.rivers > rivers.sql
-- shp2pgsql -s 3785 -d -D -I -N skip lands.shp map.lands > lands.sql
--
-- We use 'planet' database name (any valid name is suitable):
-- createdb --host=localhost --port=5432 --username=postgres \
--    --echo --template=template0 "planet" "whole planet data"
--
-- 'sectors.sql' must be evaluated first.
-- psql --host=localhost --port=5432 --username=postgres --dbname=planet /
--     --echo-errors --file=sectors.sql
-- psql --host=localhost --port=5432 --username=postgres --dbname=planet /
--     --echo-errors --file=rivers.sql
-- psql --host=localhost --port=5432 --username=postgres --dbname=planet /
--     --echo-errors --file=lands.sql

-- Optional. Cast types for some GISes (mapnik) to be worked correctly:
-- ALTER TABLE map.lands ALTER COLUMN fid TYPE BigInt;
-- ALTER TABLE map.lands ALTER COLUMN aid TYPE BigInt;
-- ALTER TABLE map.rivers ALTER COLUMN fid TYPE BigInt;
-- ALTER TABLE map.rivers ALTER COLUMN lmid TYPE BigInt;
--
-- PostgreSQL configuration:
--
-- If you intend to use large geometry data like planet data
-- from forgemaps.com, then you need to make changes
-- to PostgreSQL database configuration.
-- Look into https://postgis.net/docs/performance_tips.html
-- to understand what to do.
-- 
-- ************************************************************

CREATE EXTENSION IF NOT EXISTS postgis;

-- ************************************************************
-- Permament objects
-- ************************************************************

CREATE SCHEMA IF NOT EXISTS map;

CREATE SEQUENCE IF NOT EXISTS map."fid_gen"
    INCREMENT 1
    MINVALUE 1
    MAXVALUE 9223372036854775807
    START 1
    CACHE 1;
COMMENT ON SEQUENCE map."fid_gen"
    IS 'FID (feature id) sequence generator';

-- ocean sectors table
DROP TABLE IF EXISTS map.ocean_sectors;
CREATE TABLE IF NOT EXISTS map.ocean_sectors
(
    fid           bigint NOT NULL,
    area          double precision,
    geom          geometry(POLYGON, 3785),
    center        geometry(POINT, 3785)
)
WITH (
    OIDS=FALSE
);

-- ocean sectors graph table
DROP TABLE IF EXISTS map.ocean_sectors_graph;
CREATE TABLE IF NOT EXISTS map.ocean_sectors_graph
(
    fid_l  BigInt NOT NULL,
    fid_r  BigInt NOT NULL
)
WITH (
    OIDS=FALSE
);

-- land sectors table
DROP TABLE IF EXISTS map.land_sectors;
CREATE TABLE IF NOT EXISTS map.land_sectors (
    fid          BigInt,
    sid          BigInt,
    area         Double Precision,
    geom         geometry(MULTIPOLYGON, 3785),
    aid          BigInt
)
WITH (
    OIDS=FALSE
);

CREATE UNIQUE INDEX land_sectors_idx ON map.land_sectors (fid, sid);
CREATE INDEX land_sectors_aid_idx ON map.land_sectors (aid);
CREATE INDEX land_sectors_geom_idx ON map.land_sectors USING GIST (geom);

-- land sectors graph table
DROP TABLE IF EXISTS map.land_sectors_graph;
CREATE TABLE IF NOT EXISTS map.land_sectors_graph (
    aid    BigInt,
    fid_l  BigInt NOT NULL,
    fid_r  BigInt NOT NULL
)
WITH (
    OIDS=FALSE
);

CREATE INDEX land_sectors_graph_idx ON map.land_sectors_graph (aid, fid_l);

-- ************************************************************
-- Base functions
-- ************************************************************

-- Does fast intersection by preliminary clipping bigger geometry 'geomB'
-- with the bbox of lesser geometry 'geomA'
DROP FUNCTION IF EXISTS map.fastIntersection;
CREATE OR REPLACE FUNCTION map.fastIntersection(
    geomA Geometry, geomB Geometry )
    RETURNS Geometry AS $proc$
DECLARE
BEGIN
    RETURN ST_Intersection(geomA, ST_ClipByBox2D(geomB, ST_Envelope(geomA)));
END;
$proc$ LANGUAGE plpgsql;

-- Splits any geometry into separate polygons (if they exist)
DROP FUNCTION IF EXISTS map.splitIntoPolygons;
CREATE OR REPLACE FUNCTION map.splitIntoPolygons(geom Geometry)
    RETURNS SETOF Geometry AS $proc$
DECLARE
BEGIN
    IF (geom IS NULL) OR (ST_GeometryType(geom) = 'ST_Polygon') THEN
        RETURN NEXT geom;
    ELSE 
        RETURN QUERY SELECT a.gdump
            FROM (SELECT (ST_Dump(geom)).geom AS gdump) AS a
            WHERE ST_GeometryType(a.gdump)='ST_Polygon';
    END IF;

    RETURN;
END;
$proc$ LANGUAGE plpgsql;

-- Makes random initial points (seeds) for voronoi diagram. Points appear
-- less frequently near the poles, in accordance with the decrease
-- in the length of latitude circles
DROP FUNCTION IF EXISTS map.genPoints;
CREATE OR REPLACE FUNCTION map.genPoints(env4326 Geometry, n_points Int)
    RETURNS SETOF Geometry AS $proc$
DECLARE
    xMin Int;
    xMax Int;
    yMin Int;
    yMax Int;
    x Double Precision;
    y Double Precision;
    maxCosY Double Precision;
    count Int := 0;
BEGIN
    xMin := ST_XMin(env4326);
    xMax := ST_XMax(env4326);
    yMin := ST_YMin(env4326);
    yMax := ST_YMax(env4326);

    IF yMin*yMax<=0 THEN -- area crosses or touches the equator
        maxCosY := 1;
    ELSEIF abs(yMin)<abs(yMax) THEN -- top hemisphere
        maxCosY := cosd(yMin);
    ELSE
        maxCosY := cosd(yMax); -- bottom hemisphere
    END IF;

    LOOP
        y := (yMax-yMin) * random() + yMin;
        CONTINUE WHEN random()>=cosd(y)/maxCosY; -- discard excess points
        x := (xMax-xMin) * random() + xMin;
        count := count + 1;
        RETURN NEXT ST_SetSRID(ST_MakePoint(x, y), 4326);
        IF count >= n_points THEN
            EXIT;
        END IF;
    END LOOP;
    
    RETURN;
END;
$proc$ LANGUAGE plpgsql;

-- Construct MultiPoint geometry for voronoi diagram seeds
DROP FUNCTION IF EXISTS map.genMultiPoint;
CREATE OR REPLACE FUNCTION map.genMultiPoint(
    env4326 Geometry, avg_vp_area Double Precision, is_whole_planet Boolean )
    RETURNS Geometry AS $proc$
DECLARE
    area_unit Double Precision := 1000000;
    planet4326 Geometry := ST_GeomFromText(
        'POLYGON((-180 -89.9999, 180 -89.9999, 180 89.9999, -180 89.9999, -180 -89.9999))',
        4326 );
    points Geometry;
    n_points Int;
BEGIN
    IF is_whole_planet THEN
        n_points := ceiling(510072000.0 * area_unit / avg_vp_area);
        SELECT INTO points (
            SELECT ST_Collect(p.the_geom)
            FROM (SELECT map.genPoints(planet4326, n_points) AS the_geom) AS p
            WHERE ST_ContainsProperly(env4326, p.the_geom)
        );
    ELSE
        -- XXX: ST_Area doesn't work properly for very large areas
        n_points := ceiling(ST_Area(Geography(env4326), False) / avg_vp_area);
        SELECT INTO points (
            SELECT ST_Collect(p.the_geom)
            FROM (SELECT map.genPoints(env4326, n_points) AS the_geom) AS p
        );
    END IF;

    IF (points IS NULL) OR (ST_NPoints(points)<=1) THEN
        -- if zero or one point then don't make result geometry
        points := NULL;
        RAISE DEBUG 'seeds points: %', 0;
    ELSE 
        points := ST_Transform(points, 3785);
        RAISE DEBUG 'seeds points: %', ST_NPoints(points);
    END IF;

    RETURN points;
END;
$proc$ LANGUAGE plpgsql;

-- Makes voronoi diagram in 'env4326' geometry with 'points' seeds
-- and points placement 'tolerance'
DROP FUNCTION IF EXISTS map.makeVoronoiPolygons;
CREATE OR REPLACE FUNCTION map.makeVoronoiPolygons(
    env4326 Geometry, points Geometry, tolerance Float)
    RETURNS Geometry AS $proc$
DECLARE
    vorColl Geometry;
BEGIN
    IF points IS NOT NULL THEN
        SELECT INTO vorColl (
            SELECT ST_VoronoiPolygons(
                points, tolerance, ST_Transform(env4326, 3785) )
        );
        IF (vorColl IS NULL) OR (ST_IsEmpty(vorColl)) THEN 
            RAISE EXCEPTION '*** empty voronoi diagram';
        END IF;
    ELSE 
        -- if too few points then make only one polygon
        vorColl := ST_Transform(env4326, 3785);
    END IF;

    RETURN vorColl;
END;
$proc$ LANGUAGE plpgsql;

-- ************************************************************
-- Ocean sectors
-- ************************************************************

-- Makes ocean voronoi polygons in 'worldBB3785' area with average
-- area 'avg_area'. If 'do_relaxation' is true then makes Lloyd's relaxation
-- until stabilization over poligon minimum area is obtained
DROP FUNCTION IF EXISTS map.makeOceanVoronoiPolygons;
CREATE OR REPLACE FUNCTION map.makeOceanVoronoiPolygons(
    worldBB3785 Geometry, avg_area Double Precision, is_whole_planet Boolean,
    relax_level Int )
    RETURNS Void AS $proc$
DECLARE
    env4326 Geometry;
    tolerance Double Precision;
    vorColl Geometry;   -- voronoi polygons as GeometryCollection
    points Geometry;
    --min_area Double Precision := 0;
    --pred_min_area Double Precision := 0;
    vp_geom Geometry;
    c_relax_level Int := 0;
    curs_vps CURSOR FOR
        SELECT map.splitIntoPolygons(
            ST_Intersection(ST_SetSRID(a.gdump, 3785), worldBB3785)) as geom
        FROM (SELECT (ST_Dump(vorColl)).geom AS gdump) AS a;
    curs_lands CURSOR (geom0 Geometry) FOR
        SELECT fid, aid, area, geom FROM map.lands
        WHERE inlake=False AND geom && geom0
        ORDER BY area DESC;
BEGIN
    -- tolerance of points placement (shoud be function of 'avg_area')
    tolerance := Sqrt(avg_area)/5;

    -- temporary table for voronoi polygons
    DROP TABLE IF EXISTS _vps;
    CREATE TEMP TABLE _vps (
        fid           BigInt PRIMARY KEY DEFAULT nextval('map."fid_gen"'),
        geom          geometry(POLYGON, 3785)
    );

    -- initial MultiPoints for the Lloyd's algorithm
    env4326 := ST_Envelope(ST_Transform(worldBB3785, 4326));
    points := map.genMultiPoint(env4326, avg_area, is_whole_planet);

    -- Lloyd's algorithm iterations
    <<relaxation>>
    LOOP
        SELECT map.makeVoronoiPolygons(env4326, points, tolerance)
            INTO vorColl;
        EXIT relaxation WHEN c_relax_level >= relax_level;
        RAISE DEBUG 'relaxation level: %', c_relax_level;
        SELECT ST_Collect(ST_Centroid(
            ST_Intersection(a.gdump, worldBB3785) ))
            INTO points     -- seeds for next iteration
        FROM (SELECT (ST_Dump(vorColl)).geom AS gdump) AS a;
        c_relax_level := c_relax_level + 1;
    END LOOP relaxation;

    RAISE LOG '(o) voronoi polygons: %', ST_NumGeometries(vorColl);
    RAISE LOG '(o) substracting land areas..';

    -- substract land areas from voronoi polygons making sectors
    <<curs_vps>>
    FOR vp IN curs_vps LOOP
        vp_geom = vp.geom;
        <<curs_lands>>
        FOR land IN curs_lands (vp_geom) LOOP
            vp_geom := ST_Difference(vp_geom, land.geom);
        END LOOP curs_lands;
        INSERT INTO _vps (geom)
        SELECT map.splitIntoPolygons(vp_geom);
    END LOOP curs_vps;
END;
$proc$ LANGUAGE plpgsql;

-- Calcs effective minimum area. Make the requirement for a minimum area
-- larger for sectors near the poles. It does not compensate
-- the shrinking of hi-latitude sectors in the process of relaxation,
-- but somewhat helps prevent their ugly merging.
-- (Note. Second possible way (not inplemented) is to make initially seeds
-- for voronoi diagram even more sparsely towards the poles,
-- but the dependency of required rarefaction on the number
-- of relaxation steps is not yet clear)
DROP FUNCTION IF EXISTS map.effArea(Double Precision, Geometry);
CREATE OR REPLACE FUNCTION map.effArea(
    area Double Precision, center Geometry )
    RETURNS Double Precision AS $proc$
DECLARE
    y Double Precision;
BEGIN
    IF ST_GeometryType(center)<>'ST_Point' THEN
        RAISE EXCEPTION 'effArea: invalid type of center point: %:',
            ST_GeometryType(center);
    END IF;
    y := ST_Y(ST_Transform(center, 4326));
    RETURN cosd(y/2) * area;
END;
$proc$ LANGUAGE plpgsql;

-- Merging small sectors using 'merging_method' (currently, '1' or '2')
-- '1': merging along the longest border
-- '2': merging to the polygon with the nearest centroid
DROP FUNCTION IF EXISTS map.mergeOceanSmallSectors;
CREATE OR REPLACE FUNCTION map.mergeOceanSmallSectors(
    min_area Double Precision, merging_method Int )
    RETURNS BigInt AS $proc$
DECLARE
    weight Double Precision;
    d Double Precision;
    l Double Precision;
    count_merged BigInt := 0;
    sel_fid BigInt;
    sel_geom Geometry;
    new_geom Geometry;
    vp Record;
    curs_vps CURSOR FOR
        SELECT fid FROM _vps
        WHERE geom IS NOT NULL;
    curs_nbs CURSOR (fid0 BigInt, geom0 Geometry) FOR
        SELECT fid, geom, center FROM map.ocean_sectors
        WHERE geom && geom0 AND fid<>fid0 AND _ST_Intersects(geom, geom0);
BEGIN
    -- prepare initial sectors
    DROP INDEX IF EXISTS map.ocean_sectors_idx;
    DROP INDEX IF EXISTS map.ocean_sectors_geom_idx;
    TRUNCATE map.ocean_sectors;
    INSERT INTO map.ocean_sectors
        SELECT fid,
           ST_Area(Geography(ST_Transform(geom, 4326)), True),
           geom,
           ST_centroid(geom)
        FROM _vps
        WHERE geom IS NOT NULL;
    CREATE UNIQUE INDEX ocean_sectors_idx ON map.ocean_sectors (fid);
    CREATE INDEX ocean_sectors_geom_idx ON map.ocean_sectors USING GIST (geom);

    -- merge vp with a neighbour defined by 'merging_method'
    <<curs_vps_loop>>
    FOR _vp IN curs_vps LOOP
        SELECT * INTO vp FROM map.ocean_sectors WHERE fid = _vp.fid;
        IF vp IS NOT NULL AND vp.area < map.effArea(min_area, vp.center) THEN
            weight = NULL;
            sel_fid = NULL;
            <<curs_nbs_loop>>
            FOR nb IN curs_nbs (vp.fid, vp.geom) LOOP
                IF merging_method=2 THEN
                    d = ST_Distance(vp.center, nb.center);
                    IF (weight IS NULL) OR (d<weight) THEN
                        weight = d;
                        sel_fid = nb.fid;
                    END IF;
                ELSIF merging_method=1 THEN
                    l = ST_Length(ST_Intersection(vp.geom, nb.geom));
                    IF (weight IS NULL) OR (l>weight) THEN
                        weight = l;
                        sel_fid = nb.fid;
                    END IF;
                ELSE
                    RAISE EXCEPTION 'Unknown merging method: %', merging_method;
                END IF;
            END LOOP curs_nbs_loop;
            IF sel_fid IS NOT NULL THEN
                SELECT geom INTO sel_geom FROM map.ocean_sectors
                WHERE fid=sel_fid;
                new_geom = ST_Union(sel_geom, ST_Buffer(vp.geom, 1));
                IF ST_GeometryType(new_geom) = 'ST_Polygon' THEN
                    DELETE FROM map.ocean_sectors WHERE fid=vp.fid;
                    UPDATE map.ocean_sectors SET
                        geom=new_geom,
                        area=ST_Area(Geography(
                            ST_Transform(new_geom, 4326) ), True),
                        center=ST_centroid(new_geom)
                    WHERE fid=sel_fid;
                    count_merged = count_merged + 1;
                ELSE
                    -- XXX: currently, we leave polygons as is
                    RAISE WARNING 'Strange geometry for % (%)',
                        vp.fid, ST_GeometryType(new_geom);
                END IF;
            ELSE
                RAISE DEBUG 'Cannot find neighbour for %', vp.fid;
            END IF;
        END IF;
    END LOOP curs_vps_loop;

    --DROP TABLE _vps;
    RETURN count_merged;
END;
$proc$ LANGUAGE plpgsql;

-- Builds graph of neighbouring ocean sectors
DROP FUNCTION IF EXISTS map.buildOceanSectorsGraph;
CREATE OR REPLACE FUNCTION map.buildOceanSectorsGraph()
    RETURNS Integer AS $proc$
DECLARE
    count Int := 0;
    curs_sectors CURSOR FOR
        SELECT fid, geom FROM map.ocean_sectors;
    curs_nb CURSOR (fid0 BigInt, geom0 Geometry) FOR
        SELECT * FROM map.ocean_sectors
        WHERE geom && geom0 AND fid<>fid0 AND _ST_Intersects(geom, geom0);
BEGIN
    DROP INDEX IF EXISTS map.ocean_sectors_graph_idx;
    TRUNCATE map.ocean_sectors_graph;

    <<curs_sectors_loop>>
    FOR sector IN curs_sectors LOOP
        <<curs_nb_loop>>
        FOR nb IN curs_nb (sector.fid, sector.geom) LOOP
            INSERT INTO map.ocean_sectors_graph VALUES (sector.fid, nb.fid);
            count = count + 1;
        END LOOP curs_nb_loop;
    END LOOP curs_sectors_loop;

    CREATE INDEX ocean_sectors_graph_idx ON map.ocean_sectors_graph (fid_l);

    RETURN count;
END;
$proc$ LANGUAGE plpgsql;

-- Main function for creating ocean sectors. Args:
-- 'world': world boundig polygon in any (not empty) projection
-- 'avg_vp_areaKM': avarage voronoi polygon area in km^2
-- 'merging_ratio': share of 'avg_vp_areaKM', such that if the sector area
--     is smaller than it, then it will be merged with another sector
-- 'merging_method': method of sectors merging (sea 'mergeOceanSmallSectors'
--     description)
DROP FUNCTION IF EXISTS map.makeOceanSectors;
CREATE OR REPLACE FUNCTION map.makeOceanSectors(
    world Geometry, avg_vp_areaKM Double Precision,
    merging_ratio Double Precision, merging_method Int
    ) RETURNS Void AS $proc$
DECLARE
    area_unit Double Precision := 1000000;
    world3785 Geometry;
    n_vps Int;
    n_sectors Int;
    graph_size Int;
BEGIN
    IF (ST_GeometryType(world)<>'ST_Polygon')
      OR (ST_IsValid(world) IS NOT TRUE) THEN
        RAISE EXCEPTION '*** incorrect world: type=%, SRID=%, valid=%',
            ST_GeometryType(world), ST_SRID(world), ST_IsValid(world);
    END IF;    
    
    world3785 = ST_Transform(world, 3785);

    RAISE NOTICE '(o) making voronoi polygons..';
    PERFORM map.makeOceanVoronoiPolygons(
        world3785, avg_vp_areaKM*area_unit, True, 3 );

    SELECT count(*) INTO n_vps FROM _vps;
    RAISE NOTICE '(o) raw polygons: %', n_vps;

    RAISE NOTICE '(o) merging small sectors..';
    SELECT map.mergeOceanSmallSectors(
        merging_ratio*avg_vp_areaKM*area_unit, merging_method
        ) INTO n_sectors;
    RAISE NOTICE '(o) sectors merged: %', n_sectors;

    RAISE NOTICE '(o) building sectors graph..';
    SELECT map.buildOceanSectorsGraph() INTO graph_size;
    RAISE NOTICE '(o) sectors graph size: %', graph_size;
END;
$proc$ LANGUAGE plpgsql;

-- ************************************************************
-- Land sectors
-- ************************************************************

-- Makes lands voronoi polygons for lands with defined 'aid' and
-- average 'avg_vp_area' area. Orders voronoi polygons randomly
DROP FUNCTION IF EXISTS map.makeLandVoronoiPolygons;
CREATE OR REPLACE FUNCTION map.makeLandVoronoiPolygons(
    aid BigInt, avg_vp_area Double Precision )
    RETURNS Int AS $proc$
DECLARE
    land_part_count Int := 0;
    vorColl Geometry;
    env4326 Geometry;
    points Geometry;
    tolerance Double Precision;
    curs_lands CURSOR FOR
        SELECT ( CASE WHEN (ST_NumGeometries(geom)=1)
                    THEN ST_GeometryN(geom,1)
                    ELSE geom
                 END ) AS geom, fid FROM map.lands AS l
        WHERE l.aid=makeLandVoronoiPolygons.aid AND inlake=False;
BEGIN
    -- tolerance of points placement (shoud be function of 'avg_area'
    -- and maybe 'aid')
    tolerance := Sqrt(avg_vp_area)/20;

    -- temporary table for initially generated voronoi polygons
    DROP TABLE IF EXISTS _raw_vps;
    CREATE TEMP TABLE _raw_vps (
        geom          geometry(POLYGON, 3785)
    );

    -- make voronoi polygons for each land with 'aid'
    <<curs_lands_loop>>
    FOR land IN curs_lands LOOP
        land_part_count := land_part_count + 1;
        env4326 := ST_Envelope(ST_Transform(land.geom, 4326));
        points := map.genMultiPoint(env4326, avg_vp_area, False);
        SELECT map.makeVoronoiPolygons(env4326, points, tolerance)
            INTO vorColl;
        RAISE LOG '(%) % polygons for land %',
                  aid, ST_NumGeometries(vorColl), land.fid;
        INSERT INTO _raw_vps (geom)
            SELECT map.splitIntoPolygons(map.fastIntersection(
                ST_SetSRID(a.gdump, 3785), land.geom)) as geom
            FROM (SELECT (ST_Dump(vorColl)).geom AS gdump) AS a;
    END LOOP curs_lands_loop;

    -- temporary table for randomly ordered voronoi polygons
    DROP TABLE IF EXISTS _vps;
    CREATE TEMP TABLE _vps (
        fid           BigInt PRIMARY KEY DEFAULT nextval('map."fid_gen"'),
        geom          geometry(POLYGON, 3785),
        sid           Int,
        step          Int
    );

    -- order randomly voronoi poligons in '_vps' table
    INSERT INTO _vps (geom)
        SELECT * FROM _raw_vps ORDER BY random();
    DROP TABLE _raw_vps;
    CREATE INDEX _vps_idx ON _vps (fid);
    CREATE INDEX _vps_step_idx ON _vps (step);
    CREATE INDEX _vps_geom_idx ON _vps USING GIST (geom);

    RETURN land_part_count;
END;
$proc$ LANGUAGE plpgsql;

-- Makes the graph of neighbouring voronoi polygons
DROP FUNCTION IF EXISTS map.buildVPGraph;
CREATE OR REPLACE FUNCTION map.buildVPGraph()
    RETURNS Integer AS $proc$
DECLARE
    count Int := 0;
    curs_vps CURSOR FOR
        SELECT fid, geom FROM _vps;
    curs_nb CURSOR (fid0 BigInt, geom0 Geometry) FOR
        SELECT * FROM _vps
        WHERE geom && geom0 AND fid<>fid0 AND ST_Intersects(geom, geom0);
BEGIN
    DROP TABLE IF EXISTS _vps_graph;
    CREATE TEMP TABLE _vps_graph (
        fid_l  BigInt NOT NULL,
        fid_r  BigInt NOT NULL
    );

    <<curs_vps_loop>>
    FOR vp IN curs_vps LOOP
        <<curs_nb_loop>>
        FOR nb IN curs_nb (vp.fid, vp.geom) LOOP
            INSERT INTO _vps_graph VALUES (vp.fid, nb.fid);
            count := count + 1;
        END LOOP curs_nb_loop;
    END LOOP curs_vps_loop;

    CREATE INDEX vps_graph_idx ON _vps_graph (fid_l);
    RETURN count;
END;
$proc$ LANGUAGE plpgsql;

-- Assembles initial sectors by gathering neighbouring voronoi polygons
-- 'n_sectors': requested number of 'big' sectors
DROP FUNCTION IF EXISTS map.assembleRawSectors;
CREATE OR REPLACE FUNCTION map.assembleRawSectors(aid BigInt, n_sectors Int)
    RETURNS Int AS $proc$
DECLARE
    merge_sid Int;
    count Int := 0;
    l_count Int;
    n_vps Int;
    step_offset Int;
    c_step Int;
    max_steps Int;
    c_min_size Int;
    c_n_sectors Int;
    vp Record;
    curs_vps CURSOR FOR
        SELECT fid, geom FROM _vps WHERE step IS NULL;
    curs_nb CURSOR (fid0 BigInt, step_expect Int) FOR
        SELECT * FROM _vps
        INNER JOIN _raw_sectors ON _vps.sid=_raw_sectors.sid
        WHERE step<step_expect
          AND (_vps.fid IN (SELECT fid_r FROM _vps_graph WHERE fid_l=fid0));
BEGIN
    -- clean '_vps' table (just in case)
    UPDATE _vps SET sid=NULL, step=NULL;

    -- check function argument
    SELECT count(*) INTO n_vps FROM _vps;
    IF n_sectors>n_vps THEN
        RAISE EXCEPTION '*** assempleRawSectors: too many sectors (% from %)',
            n_sectors, n_vps;
    END IF;

    -- temporary table for initial sectors
    DROP TABLE IF EXISTS _raw_sectors;
    CREATE TABLE _raw_sectors (
        fid         BigInt,
        sid         BigInt,
        geom        geometry(MULTIPOLYGON, 3785),
        size        Int
    );

    -- make seed sectors
    c_step := 0;
    WITH copied_vps AS (
        (
        SELECT DISTINCT fid, geom FROM _vps
        INNER JOIN _vps_graph ON fid=fid_l LIMIT n_sectors
        )
        UNION ALL
        (
        SELECT fid, geom FROM _vps
        LEFT JOIN _vps_graph ON fid=fid_l
        WHERE fid_l IS NULL
        )
    ), inserted_vps AS (
        INSERT INTO _raw_sectors (fid, sid, geom, size)
        SELECT fid, (row_number() OVER ()) AS sid, ST_Multi(geom), 1
        FROM copied_vps AS a
        RETURNING fid, sid
    )
    UPDATE _vps SET sid=a.sid, step=0 FROM inserted_vps AS a
    WHERE _vps.fid = a.fid;

    SELECT count(*) INTO c_n_sectors FROM _raw_sectors;

    -- mark neighboring polygons with sector 'sid' on each step
    l_count := count;
    <<assemble_loop>>
    LOOP
        c_step := c_step+1;
        RAISE LOG '(%) merging step: % (so far: %)', aid, c_step, count;
        <<curs_vps_loop>>
        FOR vp IN curs_vps LOOP
            merge_sid := NULL;
            c_min_size := NULL;
            <<curs_nb_loop>>
            FOR nb IN curs_nb (vp.fid, c_step) LOOP
                IF (c_min_size IS NULL) OR (nb.size<c_min_size) THEN
                    merge_sid := nb.sid;
                    c_min_size := nb.size;
                END IF;
            END LOOP curs_nb_loop;
            IF merge_sid IS NOT NULL THEN
                UPDATE _vps SET
                    sid=merge_sid,
                    step=c_step
                WHERE CURRENT OF curs_vps;
                count := count+1;
            END IF;
        END LOOP curs_vps_loop;
        IF count = l_count THEN
            SELECT fid, geom INTO vp FROM _vps WHERE step IS NULL LIMIT 1;
            IF vp IS NOT NULL THEN
                RAISE LOG '(%) yet another sector', aid;
                c_n_sectors = c_n_sectors+1;
                INSERT INTO _raw_sectors (fid, sid, geom, size)
                VALUES (vp.fid, c_n_sectors, ST_Multi(vp.geom), 1);
                UPDATE _vps SET
                    sid=c_n_sectors,
                    step=0
                WHERE _vps.fid = vp.fid;
                c_step := 0;
            ELSE
                EXIT assemble_loop;
            END IF;
        END IF;
        IF c_step >= n_vps+1 THEN
            RAISE EXCEPTION '*** assembleRawSectors: infinit loop';
        END IF;
        l_count := count;
    END LOOP assemble_loop;

    -- union marked polygons into sectors
    RAISE LOG '(%) making unions..', aid;
    WITH assempled_vps AS (
        SELECT sid AS u_sid, count(sid) AS u_size, ST_Union(geom) AS u_geom
        FROM _vps
        GROUP BY sid
    )
    UPDATE _raw_sectors SET geom=ST_Multi(u_geom), size=u_size
    FROM assempled_vps
    WHERE sid=u_sid;

    DROP TABLE _vps;
    DROP TABLE _vps_graph;
    CREATE INDEX _raw_sectors_geom ON _raw_sectors USING GIST (geom);

    RETURN c_n_sectors;
END;
$proc$ LANGUAGE plpgsql;

-- Splits 'obj' geometry with 'blade' geometry if any of
-- the resulting geometries have area less then 'max_cut_area'
DROP FUNCTION IF EXISTS map.conditionalSplit;
CREATE OR REPLACE FUNCTION map.conditionalSplit(
    obj Geometry, blade Geometry, max_cut_area Double Precision )
    RETURNS SETOF Geometry AS $proc$
DECLARE
    part Record;
    splited_geom Geometry;
    is_do_nothing Boolean;
BEGIN
    splited_geom := ST_Split(obj, blade);
    SELECT every(
        ST_Area(Geography(ST_Transform(a.the_geom, 4326)), True)>max_cut_area )
        INTO is_do_nothing
    FROM (SELECT (ST_Dump(splited_geom)).geom AS the_geom) AS a;

    IF is_do_nothing IS TRUE THEN
        RETURN NEXT obj;
    ELSE
        FOR part IN SELECT (ST_Dump(splited_geom)).geom AS the_geom LOOP
            RETURN NEXT part.the_geom;
        END LOOP;
    END IF;

    RETURN;
END;
$proc$ LANGUAGE plpgsql;

-- Cuts sectors into pieces with rivers. Only rivers with
-- final streamflow > min_streamflow will be used
-- 'max_cut_area': one of the resulting geometries shoud have area less then it
DROP FUNCTION IF EXISTS map.cutSectors;
CREATE OR REPLACE FUNCTION map.cutSectors(
    aid BigInt, max_cut_area Double Precision, min_streamflow Int)
    RETURNS Void AS $proc$
DECLARE
    split_coll Geometry;
    curs_s CURSOR FOR
        SELECT geom, sid FROM _raw_sectors;
    curs_r CURSOR (geom0 Geometry) FOR
        SELECT a.fid, a.gdump AS geom FROM
          ( SELECT fid, (ST_Dump(ST_Force2D(geom))).geom AS gdump
            FROM map.rivers
            WHERE lmid=cutSectors.aid
              AND streamflow>=min_streamflow          -- ??
              AND geom0 && ST_Force2D(geom)
            ORDER BY streamflow ) AS a;
BEGIN
    -- temporary table for sector cuts
    DROP TABLE IF EXISTS _sector_cuts;
    CREATE TEMP TABLE _sector_cuts (
        fid          BigInt PRIMARY KEY DEFAULT nextval('map."fid_gen"'),
        sid          BigInt,
        area         double precision,
        geom         geometry(POLYGON, 3785)
    );

    -- make cuts
    <<curs_s_loop>>
    FOR sector IN curs_s LOOP
        RAISE LOG '(%) sector: %', aid, sector.sid;
        split_coll := ST_Multi(sector.geom);
        <<curs_r_loop>>
        FOR r IN curs_r (sector.geom) LOOP
            SELECT ST_Collect(b.d2) INTO split_coll
            FROM (SELECT map.conditionalSplit(a.d1, r.geom, max_cut_area) AS d2
                  FROM (SELECT (ST_Dump(split_coll)).geom AS d1) AS a) AS b;
        END LOOP curs_r_loop;
        INSERT INTO _sector_cuts (sid, geom, area)
            SELECT sector.sid,
                   a.gdump,
                   ST_Area(Geography(ST_Transform(a.gdump, 4326)), True)
                FROM (SELECT (ST_Dump(split_coll)).geom AS gdump) AS a;
    END LOOP curs_s_loop;

    DROP TABLE _raw_sectors;
    CREATE INDEX _sector_cuts_idx ON _sector_cuts (sid, area);
    CREATE INDEX _sector_cuts_geom_idx ON _sector_cuts USING GIST (geom);
END;
$proc$ LANGUAGE plpgsql;

-- Calculates 'weight' of sectors graph edge
DROP FUNCTION IF EXISTS map.getEdgeWeight;
CREATE OR REPLACE FUNCTION map.getEdgeWeight(geomA Geometry, geomB Geometry)
    RETURNS Double Precision AS $proc$
DECLARE
BEGIN
    RETURN ST_Area(ST_Intersection(geomA, geomB));
END;
$proc$ LANGUAGE plpgsql;

-- Builds sector cuts graph
DROP FUNCTION IF EXISTS map.buildSectorCutsGraph();
CREATE OR REPLACE FUNCTION map.buildSectorCutsGraph()
    RETURNS Integer AS $proc$
DECLARE
    cut_env Geometry;
    count Int := 0;
    curs_cuts CURSOR FOR
        SELECT fid, geom FROM _sector_cuts;
    curs_nb CURSOR (fid0 BigInt, geom0 Geometry) FOR
        SELECT * FROM _sector_cuts
        WHERE geom && geom0 AND fid<>fid0 AND _ST_Intersects(geom, geom0);
BEGIN
    -- temporary table for weighted cuts graph
    DROP TABLE IF EXISTS _sector_cuts_graph;
    CREATE TEMP TABLE _sector_cuts_graph (
        id     BigInt PRIMARY KEY DEFAULT nextval('map."fid_gen"'),
        fid_l  BigInt NOT NULL,
        fid_r  BigInt NOT NULL,
        weight Double Precision NOT NULL
    );

    <<curs_cuts_loop>>
    FOR cut IN curs_cuts LOOP
        <<curs_nb_loop>>
        FOR nb IN curs_nb (cut.fid, cut.geom) LOOP
            cut_env = ST_Envelope(cut.geom);
            -- we use bboxes instead of geometries for better performance
            INSERT INTO _sector_cuts_graph (fid_l, fid_r, weight) VALUES (
                cut.fid,
                nb.fid,
                map.getEdgeWeight(cut_env, nb.geom) );
            count := count + 1;
        END LOOP curs_nb_loop;
    END LOOP curs_cuts_loop;

    CREATE INDEX _sector_cuts_graph_l_idx ON _sector_cuts_graph (fid_l, weight);
    RETURN count;
END;
$proc$ LANGUAGE plpgsql;

-- Merges sector cuts
DROP FUNCTION IF EXISTS map.mergeSectorCuts;
CREATE OR REPLACE FUNCTION map.mergeSectorCuts(aid BigInt)
    RETURNS Void AS $proc$
DECLARE
    is_last_run Boolean := False;
    skip_count Int := 0;
    c_skip_count Int;
    sel_fid BigInt;
    sel_weight Double Precision;
    curs_cuts CURSOR FOR
        SELECT fid, geom, area, sid FROM _sector_cuts;
    curs_edges CURSOR (fid0 BigInt) FOR
        SELECT * FROM _sector_cuts_graph sg 
        INNER JOIN _raw_land_sectors ls ON fid_r=ls.fid
        WHERE fid_l=fid0 AND weight>0;
BEGIN
    -- temporary table for land sectors
    DROP TABLE IF EXISTS _raw_land_sectors;
    CREATE TEMP TABLE _raw_land_sectors (
        fid          BigInt,
        sid          BigInt,
        area         Double Precision,
        geom         geometry(MULTIPOLYGON, 3785)
    );

    -- move 'big' sector cuts to sectors
    WITH big_sectors AS (
        DELETE FROM _sector_cuts
        WHERE fid IN (
            SELECT DISTINCT ON (sid) fid FROM _sector_cuts
            ORDER BY sid, area DESC )
        RETURNING *
    )
    INSERT INTO _raw_land_sectors
    SELECT fid, sid, area, ST_Multi(geom) FROM big_sectors;

    -- repeat covering missed islands
    <<runs_loop>>
    LOOP
        -- merge cuts with the sector on the most 'weighted' graph edge
        c_skip_count := 0;
        <<curs_cuts_loop>>
        FOR sector_cut IN curs_cuts LOOP
            sel_fid = NULL;
            sel_weight = NULL;
            <<curs_edges_loop>>
            FOR edges IN curs_edges (sector_cut.fid) LOOP
                IF (edges.fid IS NOT NULL)
                AND ((sel_weight IS NULL) OR (sel_weight<edges.weight)) THEN
                    IF sector_cut.sid != edges.sid OR is_last_run IS TRUE THEN
                        sel_weight = edges.weight;
                        sel_fid = edges.fid;
                    END IF;
                END IF;
            END LOOP curs_edges_loop;
            IF sel_fid IS NOT NULL THEN
                WITH change_graph_l AS (
                    UPDATE _sector_cuts_graph SET
                        fid_l = sel_fid
                    WHERE fid_l=sector_cut.fid
                ), change_graph_r AS (
                    UPDATE _sector_cuts_graph SET
                        fid_r = sel_fid
                    WHERE fid_r=sector_cut.fid
                )
                UPDATE _raw_land_sectors SET
                    geom = ST_Multi(
                        ST_Union(geom, ST_Buffer(sector_cut.geom, 1)) ),
                    area = area + sector_cut.area
                WHERE fid=sel_fid;
                DELETE FROM _sector_cuts
                WHERE CURRENT OF curs_cuts;
            ELSE
                c_skip_count = c_skip_count+1;
            END IF;
        END LOOP curs_cuts_loop;

        -- check if all done
        IF (c_skip_count=0)
          OR (is_last_run IS TRUE AND c_skip_count=skip_count) THEN
            EXIT runs_loop;
        END IF;
        -- check if last run required (merging with the sector with same 'sid')
        IF is_last_run IS FALSE AND c_skip_count=skip_count THEN
            is_last_run := True;
            CONTINUE runs_loop;
        END IF;

        skip_count := c_skip_count;
    END LOOP runs_loop;
END;
$proc$ LANGUAGE plpgsql;

-- merges island with mainland sectors (if the condition is satisfied)
DROP FUNCTION IF EXISTS map.mergeIslands;
CREATE OR REPLACE FUNCTION map.mergeIslands(
    aid BigInt, pref_min_island_area Double Precision)
    RETURNS BigInt AS $proc$
DECLARE
    count_small_island Int := 0;
    intersection Record;
    container Geometry;
    r Double Precision;
    sel_in_a Double Precision;
    sel_fid BigInt;
    curs_small_islands CURSOR FOR
        SELECT fid, geom, area, sid FROM _raw_land_sectors;
BEGIN
    -- clean result table
    DELETE FROM map.land_sectors AS ls
    WHERE ls.aid=mergeIslands.aid;

    -- move nonisland sectors to result table
    WITH nonisland_rows AS (
        DELETE FROM _raw_land_sectors
        WHERE fid IN (
            SELECT DISTINCT fid FROM _raw_land_sectors AS rls
            INNER JOIN _sector_cuts_graph AS sg ON rls.fid=sg.fid_l
        )
        RETURNING *
    )
    INSERT INTO map.land_sectors
    SELECT fid, sid, area, ST_Multi(geom), aid FROM nonisland_rows;

    -- move big island to result table
    WITH big_island_rows AS (
        DELETE FROM _raw_land_sectors
        WHERE area>=pref_min_island_area
        RETURNING *
    )
    INSERT INTO map.land_sectors
    SELECT fid, sid, area, ST_Multi(geom), aid FROM big_island_rows;

    -- process small islands
    <<curs_small_islands_loop>>
    FOR small_island IN curs_small_islands LOOP
        -- try to merge small islands with sectors depending on distance
        r = 2*SQRT((pref_min_island_area-small_island.area)/PI());
        container = ST_Buffer(ST_Centroid(small_island.geom), r);
        sel_in_a = NULL;
        sel_fid = NULL;
        FOR intersection IN (
                SELECT fid, ST_Area(ST_Intersection(container, geom)) AS in_a
                FROM map.land_sectors AS ls
                WHERE ls.aid=mergeIslands.aid AND container && geom ) LOOP
            IF (sel_in_a IS NULL) OR (intersection.in_a>sel_in_a) THEN
                sel_fid = intersection.fid;
                sel_in_a = intersection.in_a;
            END IF;
        END LOOP;

        IF sel_fid IS NOT NULL THEN
            RAISE LOG '(%) merging % with %', aid, small_island.fid, sel_fid;
            UPDATE map.land_sectors AS ls SET
                geom = ST_Multi(ST_Union(ls.geom, small_island.geom)),
                area = ls.area + small_island.area
            WHERE ls.fid = sel_fid;
            DELETE FROM _raw_land_sectors
            WHERE CURRENT OF curs_small_islands;
        ELSE
            RAISE LOG '(%) separate sector: %', aid, small_island.sid;
            INSERT INTO map.land_sectors
            VALUES (
                small_island.fid,
                small_island.sid,
                small_island.area,
                ST_Multi(small_island.geom),
                aid );
            DELETE FROM _raw_land_sectors
            WHERE CURRENT OF curs_small_islands;
        END IF;

        count_small_island := count_small_island+1;
    END LOOP curs_small_islands_loop;

    -- uncomment for additional polygon simplification
    -- UPDATE map.land_sectors AS ls SET
    --     geom = map.findMultiPolygon(ST_SimplifyVW(geom, 1))
    -- WHERE ls.aid = mergeIslands.aid;

    RETURN count_small_island;
END;
$proc$ LANGUAGE plpgsql;

-- Builds final sectors graph
DROP FUNCTION IF EXISTS map.makeSectorGraph;
CREATE OR REPLACE FUNCTION map.makeSectorGraph(aid BigInt)
    RETURNS Void AS $proc$
DECLARE
BEGIN
    DELETE FROM map.land_sectors_graph AS lsg
    WHERE lsg.aid=makeSectorGraph.aid;
    
    DELETE FROM _sector_cuts_graph
    WHERE id IN ( SELECT id FROM 
          (SELECT id, ROW_NUMBER() OVER(
               PARTITION BY fid_l, fid_r ORDER BY fid_l, weight DESC) AS rn
           FROM _sector_cuts_graph) AS t
           WHERE t.rn>1 );

    INSERT INTO map.land_sectors_graph
    SELECT aid, fid_l, fid_r FROM _sector_cuts_graph;

    DROP TABLE _sector_cuts_graph;
END;
$proc$ LANGUAGE plpgsql;

-- Main function for creating land sectors.
-- 'aid': aid (areaId) of landmass
-- 'avg_vp_areaKM': avarage voronoi polygon area in km^2
-- 'avg_sector_areaKM': avarage sector area in km^2
-- 'max_sector_cut_area_ratio': share of 'avg_sector_areaKM', which defines
--     maximum area which can be cut by rivers (this is not a hard limit)
-- 'pref_min_island_area_ratio': share of 'avg_sector_areaKM', such that
--     if the island area is bigger than it, the island becomes separate sector
--     automaticaly, regardless of the distance to other sectors
-- 'min_streamflow': minimum river streamflow at which the river
--     will make the sector cuts
DROP FUNCTION IF EXISTS map.makeLandSectors;
CREATE OR REPLACE FUNCTION map.makeLandSectors(
    aid BigInt, avg_vp_areaKM Double Precision,
    avg_sector_areaKM Double Precision,
    max_sector_cut_area_ratio Float,
    pref_min_island_area_ratio Float, min_streamflow Int )
    RETURNS Void AS $proc$
DECLARE
    area_unit Int := 1000000;
    n_land_parts Int;
    n_sectors Int;
    n_vps Int;
    n_small_islands Int;
    graph_size Int;
    full_area Double Precision;
BEGIN
    IF avg_vp_areaKM>avg_sector_areaKM THEN
        RAISE EXCEPTION '(%) invalid avg_sector_areaKM and avg_vp_areaKM', aid;
    END IF;
    IF max_sector_cut_area_ratio>=0.5 THEN
        RAISE EXCEPTION '(%) invalid avg_sector_cut_area_ratio', aid;
    END IF;

    SELECT SUM(ST_Area(Geography(ST_Transform(geom, 4326)), True))
        INTO full_area FROM map.lands l WHERE l.aid=makeLandSectors.aid;
    RAISE NOTICE '(%) area: %', aid, (full_area/area_unit);
    n_vps = ceiling( full_area / avg_vp_areaKM / area_unit );
    RAISE NOTICE
        '(%) expected number of voronoi polygons: %', aid, n_vps;

    RAISE NOTICE '(%) making voronoi polygons..', aid;
    SELECT map.makeLandVoronoiPolygons(aid, avg_vp_areaKM*area_unit)
        INTO n_land_parts;
    RAISE NOTICE '(%) land parts: %', aid, n_land_parts;

    SELECT count(*) INTO n_vps FROM _vps;
    RAISE NOTICE '(%) voronoi polygons: %', aid, n_vps;

    RAISE NOTICE '(%) building voronoi graph..', aid;
    SELECT map.buildVPGraph() INTO graph_size;
    RAISE NOTICE '(%) voronoi graph size: %', aid, graph_size;

    n_sectors = ceiling( full_area / avg_sector_areaKM / area_unit );
    RAISE NOTICE
        '(%) expected number of sectors: %', aid, n_sectors;

    RAISE NOTICE '(%) assemble raw sectors..', aid;
    SELECT map.assembleRawSectors(aid, n_sectors) INTO n_sectors;
    RAISE NOTICE '(%) number of raw sectors: %', aid, n_sectors;

    RAISE NOTICE '(%) cutting sectors..', aid;
    PERFORM map.cutSectors(
        aid,
        max_sector_cut_area_ratio*avg_sector_areaKM*area_unit,
        min_streamflow
        );

    RAISE NOTICE '(%) building sectors graph..', aid;
    SELECT map.buildSectorCutsGraph() INTO graph_size;
    RAISE NOTICE '(%) graph size: %', aid, graph_size;
 
    RAISE NOTICE '(%) merging sector cuts..', aid;
    PERFORM map.mergeSectorCuts(aid);

    RAISE NOTICE '(%) merging islands..', aid;
    SELECT map.mergeIslands(
        aid, pref_min_island_area_ratio*avg_sector_areaKM*area_unit )
        INTO n_small_islands;
    RAISE NOTICE '(%) small islands: %', aid, n_small_islands;

    RAISE NOTICE '(%) making sectors graph..', aid;
    PERFORM map.makeSectorGraph(aid);
    RAISE NOTICE '(%) all done.', aid;
END;
$proc$ LANGUAGE plpgsql;
