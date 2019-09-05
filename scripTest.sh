#!/bin/bash
./histComp 133 219 ../datasets/alertas.shp 126 130 ../datasets/hidrografia.shp ../datasets/test/alertahidro.shp
./histComp 133 219 ../datasets/alertas.shp 77 75 ../datasets/rodovia.shp ../datasets/test/alertrodovia.shp
./histComp 133 219 ../datasets/alertas.shp 14 13 ../datasets/municipios.shp ../datasets/test/alertmunicipio.shp
./histComp 133 219 ../datasets/alertas.shp 7 7 ../datasets/vegeta.shp ../datasets/test/alertvegeta.shp
./histComp 126 130 ../datasets/hidrografia.shp 77 75 ../datasets/rodovia.shp ../datasets/test/hidrorodovia.shp
./histComp 126 130 ../datasets/hidrografia.shp 14 13 ../datasets/municipios.shp ../datasets/test/hidromunicipio.shp
./histComp 126 130 ../datasets/hidrografia.shp 7 7 ../datasets/vegeta.shp ../datasets/test/hidrovegeta.shp
./histComp 77 75 ../datasets/rodovia.shp 14 13 ../datasets/municipios.shp ../datasets/test/rodoviamunicipio.shp
./histComp 77 75 ../datasets/rodovia.shp 7 7 ../datasets/vegeta.shp ../datasets/test/rodoviavegeta.shp
./histComp 14 13 ../datasets/municipios.shp 7 7 ../datasets/vegeta.shp ../datasets/test/municipiovegeta.shp

./histComp 233 125 ../datasets/Rivers.shp 249 118 ../datasets/Railroads.shp ../datasets/test/riosferrovias.shp
./histComp 233 125 ../datasets/Rivers.shp 214 137 ../datasets/hydroinland.shp ../datasets/test/rioshidroinland.shp
./histComp 233 125 ../datasets/Rivers.shp 181 161 ../datasets/elev-contour-l.shp ../datasets/test/rioscontorno.shp
./histComp 233 125 ../datasets/Rivers.shp 203 78 ../datasets/Crops.shp ../datasets/test/rioscultura.shp
./histComp 249 118 ../datasets/Railroads.shp 214 137 ../datasets/hydroinland.shp ../datasets/test/ferroviashydroinland.shp
./histComp 249 118 ../datasets/Railroads.shp 181 161 ../datasets/elev-contour-l.shp ../datasets/test/ferroviascontorno.shp
./histComp 249 118 ../datasets/Railroads.shp 203 78 ../datasets/Crops.shp ../datasets/test/ferroviascultura.shp
./histComp 214 137 ../datasets/hydroinland.shp 181 161 ../datasets/elev-contour-l.shp ../datasets/test/hydroinlandcontorno.shp
./histComp 214 137 ../datasets/hydroinland.shp 203 78 ../datasets/Crops.shp ../datasets/test/hydroinlandcultura.shp
./histComp 181 161 ../datasets/elev-contour-l.shp 203 78 ../datasets/Crops.shp ../datasets/test/contornocultura.shp

