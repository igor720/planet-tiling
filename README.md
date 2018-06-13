# planet-tiling
Polygons tiling of planet sphere with PostGIS

##Data preparation

The following commands suggest default PostgreSQL+PostGIS
installation. 'riversz' and 'lands' shapefiles should be fetched
from ForgedMaps.com or they must contain data with similar structure

Prepare import sql files (we use 'map' schema for planet map data):

```
shp2pgsql -s 3785 -d -D -I -t 2D -N skip riversz.shp map.rivers > rivers.sql
shp2pgsql -s 3785 -d -D -I -N skip lands.shp map.lands > lands.sql
```

We use 'planet' database name (any valid name is suitable):

```
createdb --host=localhost --port=5432 --username=postgres \
    --echo --template=template0 "planet" "whole planet data"

'sectors.sql' must be evaluated first.

```
psql --host=localhost --port=5432 --username=postgres --dbname=planet /
    --echo-errors --file=sectors.sql
psql --host=localhost --port=5432 --username=postgres --dbname=planet /
    --echo-errors --file=rivers.sql
psql --host=localhost --port=5432 --username=postgres --dbname=planet /
    --echo-errors --file=lands.sql
```

Optional. Cast types for some GISes (mapnik) to be worked correctly:

```sql
ALTER TABLE map.lands ALTER COLUMN fid TYPE BigInt;
ALTER TABLE map.lands ALTER COLUMN aid TYPE BigInt;
ALTER TABLE map.rivers ALTER COLUMN fid TYPE BigInt;
ALTER TABLE map.rivers ALTER COLUMN lmid TYPE BigInt;
```

##Synopsis

Generate ocean sectors
```sql
map.makeOceanSectors(
    world Geometry,
    avg_vp_areaKM Double Precision,
    merging_ratio Double Precision,
    merging_method Int,
    is_whole_planet Boolean
    ) RETURNS Void
```

Generate land sectors
```sql
map.makeLandSectors(
    aid BigInt,
    avg_vp_areaKM Double Precision,
    avg_sector_areaKM Double Precision,
    max_sector_cut_area_ratio Float,
    pref_min_island_area_ratio Float,
    min_streamflow Int
    ) RETURNS Void
```

##Examples

Following will generate ocean tiles in the defined polygon with average
size of 1000000 km^2, merge small sectors if thier size less then half
of average size, using 'second' merging method and taking into account
only the specified world polygon for calculating the number of voronoi polygons.

```sql
SELECT * FROM map.makeOceanSectors(
     ST_GeomFromText('POLYGON((-75 -85, 75 -85, 75 85, -75 85, -75 -85))', 4326),
     1000000, 0.5, 2, False );
```

Following will generate land tiles within specified land (aid=5) with average
voronoi poligons size of 5000 km^2, with average sector size of 40000 km^2,
with maximum area of 0.125*40000 km^2 which can be cutted by rivers,
with minimum area of 0.25*40000 km^2 intended for sector-independent islands,
using only rivers with minimum streamflow of 2 for sectors cutting.

```sql
SELECT * FROM map.makeLandSectors(5, 5000, 40000, 0.125, 0.25, 2);
```


