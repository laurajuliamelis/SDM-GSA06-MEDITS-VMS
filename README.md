[![DOI](https://zenodo.org/badge/593249638.svg)](https://zenodo.org/badge/latestdoi/593249638)

# SDM-GSA06-MEDITS-VMS: Código y datos del Proyecto SteLar-DELTA

## Información general 
Este repositorio contiene el *script* de R y los resultados del análisis espacial llevado a cabo en el marco del contrato de apoyo tecnológico firmado entre la Agencia Estatal Consejo Superior de Investigaciones Científicas (CSIC) y CORY's - Investigación y Conservación de la Biodiversidad. El informe de resultados puede enontrarse [aquí](https://digital.csic.es/handle/10261/295428).

Se ha evaluado la distribución espacial en términos de biomasa y probabilidad de ocurrencia de las tres principales especies de pequeños peces pelágicos (anchoa, sardina y alacha) en la plataforma continental ibérica dentro de la sub-área geográfica (GSA) 06, y particularmente dentro de la Zona de Especial Protección para las Aves (ZEPA) ES0000512.

Para hacer estas estimas combinamos dos bases de datos independientes (procedentes de campañas oceanográficas MEDITS y datos de capturas reportados en los diarios de pesca de buques arrastreros) con información ambiental dentro de avanzados modelos de distribución de especies (BRT-RAC).

## Contenido

En el repositorio se incluye tanto el código utilizado para llevar a cabo el análisis como los resultados obtenidos.

### [**Figuras.**](Figuras)

La carpeta contiene las figuras incluidas en el informe del proyecto: área de estudio, mapas de las variables ambientales, tendencias temporales y funciones de densidad de probabilidad de la biomasa de las tres especies, mapas de probabilidad de ocurrencia predicha para cada especie (anual -media de junio y julio- y media de los años 1998 a 2019) y mapas de biomasa predicha para cada especie (anual -media de junio y julio- y media de los años 1998 a 2019).

### [**Predicciones.**](Predicciones)

Contiene las predicciones de probabilidad de ocurrencia y de biomasa para cada especie (anchoa, sardina y alacha) para los años 1998 a 2019 (media de junio y julio), así como la media de todos los años y su desviación estándar.

### [**R.**](R)

El archivo `BRT-RAC_Hurdle-Models.R` contiene el código utilizado para ajustar modelos BRT-RAC (*Boosted Regression Trees - Residual Autocovariate*) con el método delta log-normal o *hurdle model*.

Los documentos `gadm36_ESP_0_sp.rds` y `gadm36_FRA_0_sp.rds` contienen los polígonos de España y Francia, respectivamente. Estos se han obenido mediante la función `getData()` del paquete <code>&lt;raster&gt;</code> (version 3.6-14).

### [**Shapefiles GSA-ZEPA.**](Shapefiles%20GSA-ZEPA)

Archivos en formato simple con las áreas (i.e. polígonos) de la sub-área geográfica (GSA) 06 y la Zona de Especial Protección para las Aves del Delta del Ebro (ZEPA ES0000512).


## Licencia

Este proyecto está bajo la licencia MIT - para más detalles, consulte el archivo [LICENSE.md](LICENSE.md)

