command line wrapper for WorldWind elevation API

1. export as runnable jar with 'package required libraries' option in eclipse
2. put into main WorldWind folder e.g. x:\java\worldwind
3. launch with:
java -Xmx1024m -Dsun.java2d.noddraw=true -cp x:\java\worldwind\WorldWindAPI.jar;x:\java\worldwind\worldwind.jar;x:\java\worldwind\worldwindx.jar;x:\java\worldwind\jogl-all.jar;x:\java\worldwind\gluegen-rt.jar;x:\java\worldwind\gdal.jar teleworx.worldwind.WorldWindAPI test.csv
