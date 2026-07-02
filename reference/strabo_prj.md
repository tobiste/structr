# Example of StraboSpot Project

Strabospot poject containing 347 foliation and stretching lineation
measurements from the Shebandowan Greenstone Belt and Quetico
Subprovince (Stephan et al., 2025)

## Usage

``` r
data('strabo_prj')
```

## Format

An object of class `strabo`

## References

Stephan, T., Phillips, N., Tiitto, H., Perez, A., Nwakanma, M., Creaser,
R., & Hollings, P. (2025). Going with the flow - Changes of vorticity
control gold enrichment in Archean shear zones (Shebandowan Greenstone
Belt, Superior Province, Canada). Journal of Structural Geology, 201,
105542.
[doi:10.1016/j.jsg.2025.105542](https://doi.org/10.1016/j.jsg.2025.105542)

## See also

Other datasets:
[`angelier1990`](https://tobiste.github.io/structr/reference/angelier1990.md),
[`example_lines`](https://tobiste.github.io/structr/reference/example_lines.md),
[`example_lines_df`](https://tobiste.github.io/structr/reference/example_lines_df.md),
[`example_planes`](https://tobiste.github.io/structr/reference/example_planes.md),
[`example_planes_df`](https://tobiste.github.io/structr/reference/example_planes_df.md),
[`gray_example`](https://tobiste.github.io/structr/reference/gray_example.md),
[`holst`](https://tobiste.github.io/structr/reference/holst.md),
[`hossack1968`](https://tobiste.github.io/structr/reference/hossack1968.md),
[`osmundsen2010`](https://tobiste.github.io/structr/reference/osmundsen2010.md),
[`ramsay`](https://tobiste.github.io/structr/reference/ramsay.md),
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md)

## Examples

``` r
data("strabo_prj")
head(strabo_prj)
#> $data
#>                                        id dip_direction   dip strike trend
#>                                    <char>         <num> <num>  <num> <num>
#>   1: d951938c-6bee-4ce4-bee7-6cf17693a090           332    87    242    62
#>   2: 4a0aba4f-4707-4852-9811-f8457fd226e7           337    89    247    67
#>   3: 86c34f8d-b277-4a27-8f49-9339df976f66           145    84     55   235
#>   4: ec15d95a-e23c-4ba5-8280-63019e7ce553           330    82    240    60
#>   5: ad910789-d6b2-4179-9c38-9e95149d7b56           334    81    244   246
#>  ---                                                                      
#> 319: 3ae2760c-136b-4cb4-bc37-95c852642185           325    82    235   237
#> 320: 2810a725-9eee-47d7-9c97-5b5812092240           314    73    224    44
#> 321: 88d22c99-21ad-460d-b364-a8bd6dab7c7e           158    90     68   248
#> 322: 64086f81-b490-40f7-8291-f37d8c903d41           347    88    257   258
#> 323: e9b4bca6-37d7-4a2b-8833-ddc711da0a79           320    82    230   232
#>      plunge associated planar_type quality unix_timestamp  notes
#>       <num>     <lgcl>      <char>  <char>          <num> <char>
#>   1:      8       TRUE        <NA>       5   1.665089e+12   <NA>
#>   2:      4       TRUE        <NA>       5   1.665089e+12   <NA>
#>   3:      2       TRUE        <NA>       5   1.665089e+12   <NA>
#>   4:      0       TRUE        <NA>       5   1.665089e+12   <NA>
#>   5:     14       TRUE        <NA>       5   1.665089e+12   <NA>
#>  ---                                                            
#> 319:     12       TRUE        <NA>       4   1.721928e+12   <NA>
#> 320:      1       TRUE        <NA>       5   1.721928e+12   <NA>
#> 321:     15       TRUE        <NA>       3   1.721928e+12   <NA>
#> 322:      4       TRUE        <NA>       3   1.721928e+12   <NA>
#> 323:     14       TRUE        <NA>       5   1.721930e+12   <NA>
#>      modified_timestamp          spot        spot_id feature_type
#>                   <num>        <char>         <char>       <char>
#>   1:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   2:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   3:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   4:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   5:                 NA 22-TS-Moss-02 16650869118218    foliation
#>  ---                                                             
#> 319:                 NA       24TS-16 17219277194544    foliation
#> 320:                 NA       24TS-16 17219277194544    foliation
#> 321:                 NA       24TS-16 17219277194544    foliation
#> 322:                 NA       24TS-16 17219277194544    foliation
#> 323:                 NA       24TS-17 17219297078328    foliation
#>      foliation_type         label contact_type fault_or_sz_type bedding_type
#>              <char>        <char>       <char>           <char>       <char>
#>   1:       cleavage foliation 242         <NA>             <NA>         <NA>
#>   2:       cleavage foliation 247         <NA>             <NA>         <NA>
#>   3:       cleavage  foliation 55         <NA>             <NA>         <NA>
#>   4:       cleavage foliation 240         <NA>             <NA>         <NA>
#>   5:       cleavage foliation 244         <NA>             <NA>         <NA>
#>  ---                                                                        
#> 319:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 320:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 321:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 322:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 323:           <NA>          <NA>         <NA>             <NA>         <NA>
#>      vein_fill directional_indicators vein_type other_vein_fill
#>         <char>                 <char>    <char>          <char>
#>   1:      <NA>                   <NA>      <NA>            <NA>
#>   2:      <NA>                   <NA>      <NA>            <NA>
#>   3:      <NA>                   <NA>      <NA>            <NA>
#>   4:      <NA>                   <NA>      <NA>            <NA>
#>   5:      <NA>                   <NA>      <NA>            <NA>
#>  ---                                                           
#> 319:      <NA>                   <NA>      <NA>            <NA>
#> 320:      <NA>                   <NA>      <NA>            <NA>
#> 321:      <NA>                   <NA>      <NA>            <NA>
#> 322:      <NA>                   <NA>      <NA>            <NA>
#> 323:      <NA>                   <NA>      <NA>            <NA>
#>      foliation_defined_by movement other_feature fracture_type
#>                    <char>   <char>        <char>        <char>
#>   1:                 <NA>     <NA>          <NA>          <NA>
#>   2:                 <NA>     <NA>          <NA>          <NA>
#>   3:                 <NA>     <NA>          <NA>          <NA>
#>   4:                 <NA>     <NA>          <NA>          <NA>
#>   5:                 <NA>     <NA>          <NA>          <NA>
#>  ---                                                          
#> 319:                 <NA>     <NA>          <NA>          <NA>
#> 320:                 <NA>     <NA>          <NA>          <NA>
#> 321:                 <NA>     <NA>          <NA>          <NA>
#> 322:                 <NA>     <NA>          <NA>          <NA>
#> 323:                 <NA>     <NA>          <NA>          <NA>
#>      movement_justification        linear_type linear_quality linear_notes
#>                      <char>             <char>         <char>       <char>
#>   1:                   <NA> linear_orientation              5         <NA>
#>   2:                   <NA> linear_orientation              5         <NA>
#>   3:                   <NA> linear_orientation              5         <NA>
#>   4:                   <NA> linear_orientation              5         <NA>
#>   5:                   <NA> linear_orientation              5         <NA>
#>  ---                                                                      
#> 319:                   <NA> linear_orientation              4         <NA>
#> 320:                   <NA> linear_orientation              5         <NA>
#> 321:                   <NA> linear_orientation              3         <NA>
#> 322:                   <NA> linear_orientation              3         <NA>
#> 323:                   <NA> linear_orientation              5         <NA>
#>      linear_feature_type linear_defined_by   linear_label               type
#>                   <char>            <char>         <char>             <char>
#>   1:          stretching              <NA>  stretching 62 planar_orientation
#>   2:          stretching              <NA>  stretching 67 planar_orientation
#>   3:          stretching              <NA> stretching 235 planar_orientation
#>   4:          stretching              <NA>  stretching 60 planar_orientation
#>   5:          stretching              <NA> stretching 246 planar_orientation
#>  ---                                                                        
#> 319:       mineral_align              <NA>           <NA> planar_orientation
#> 320:       mineral_align              <NA>           <NA> planar_orientation
#> 321:       mineral_align              <NA>           <NA> planar_orientation
#> 322:       mineral_align              <NA>           <NA> planar_orientation
#> 323:       mineral_align              <NA>           <NA> planar_orientation
#>      linear_unix_timestamp                            linear_id defined_by
#>                      <num>                               <char>     <char>
#>   1:          1.665089e+12 50b1dbef-5d64-4a8c-b427-7b1bb8c6053e       <NA>
#>   2:          1.665089e+12 af52e337-9c6b-4392-9f75-eb2dea09e016       <NA>
#>   3:          1.665089e+12 8189eda4-3acd-4bb8-9075-7cf5b6e3b2c1       <NA>
#>   4:          1.665089e+12 c8d4b893-6d70-4f89-b46b-ace81a01136d       <NA>
#>   5:          1.665089e+12 bf08e531-8c0e-4337-b9a3-7eccd8019486       <NA>
#>  ---                                                                      
#> 319:          1.721928e+12 84bcb839-a8fa-4a25-b28e-16b97ecbdce6       <NA>
#> 320:          1.721928e+12 9470b29b-550a-4613-9834-8667684f2aa6       <NA>
#> 321:          1.721928e+12 ceb03667-90cd-45bd-adfe-24102cbe7e42       <NA>
#> 322:          1.721928e+12 e55af915-d9a1-4d38-a710-71172ef34455       <NA>
#> 323:          1.721930e+12 9ff0f807-d54e-4dd9-b79c-57f1a1ff1116       <NA>
#>      linear_other_feature
#>                    <char>
#>   1:                 <NA>
#>   2:                 <NA>
#>   3:                 <NA>
#>   4:                 <NA>
#>   5:                 <NA>
#>  ---                     
#> 319:                 <NA>
#> 320:                 <NA>
#> 321:                 <NA>
#> 322:                 <NA>
#> 323:                 <NA>
#> 
#> $spots
#> Simple feature collection with 93 features and 22 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -90.95092 ymin: 48.44102 xmax: -90.07832 ymax: 48.67911
#> Geodetic CRS:  WGS 84
#> First 10 features:
#>            spot_id          ds_id          spot longitude latitude altitude
#> 1   16650869118218 16649847460086 22-TS-Moss-02 -90.71200 48.54260 433.6175
#> 5   16652487703646 16649847460086 22-TS-Moss-06 -90.72518 48.53076 433.4179
#> 9   16652618160235 16649847460086 22-TS-Moss-10 -90.72712 48.57296 451.8851
#> 15  16653282590279 16649847460086 22-TS-Moss-13 -90.58906 48.60464 465.7021
#> 17  16653332740820 16649847460086 22-TS-Moss-15 -90.58110 48.60599 462.3942
#> 24  16653468274092 16649847460086 22-TS-Moss-18 -90.57892 48.56819 470.4441
#> 32  16654149768459 16649847460086 22-TS-Moss-22 -90.60366 48.54289 483.4684
#> 146 16654974392162 16649847460086 22-TS-Moss-31 -90.74599 48.51911 455.9059
#> 152 16655010610405 16649847460086 22-TS-Moss-36 -90.76746 48.52015 452.1578
#> 165 16655121372131 16649847460086 22-TS-Moss-41 -90.79470 48.53087 443.9397
#>                                                                                                                                                                          note
#> 1                                                   Diorite, deformed, sc’ fabric, sinistral shear sense, scaly fabric with phacoids ranging from <1cm to several cm in scale
#> 5                                                                              Diorite and sc diorite, dextral sense of shear, along w-e shear zones\nsketch in TS field book
#> 9                                                                                                                                                          ?diorite, deformed
#> 15                                                                                                                  Sc-Diorite, sc’ shear bands, sinistral along sw-ne zones.
#> 17                                                                                                                                                                       <NA>
#> 24                                                      Greenschist, strong cleavage, ?sinistral\nFelsic metavolcanic rock, ?meta-tuffite\nKink bands with vertical fold axis
#> 32                                                                                                                                         ?Grt-Amphibolite, gneissic texture
#> 146                                                         Greenschist, cleavage, disseminated sulfides in matrix\nMassive sulfides /w Qz along veins, ?parallel to cleavage
#> 152                                                                                                                          Greenschist, ?felsic metavolcanic rock, foliated
#> 165 Greenschist, ?mafic metavolcanic rock. Foliated (cleavage). Rotated Qtz and Py-clast show sinistral shear along NE-SW. Opposite sense observed (indicating general shear)
#>                    time gps_accuracy ds_name tag:Diabase dikes tag:Diorite
#> 1   2022-10-06 20:08:31     4.603851      TS             FALSE        TRUE
#> 5   2022-10-08 17:06:10     4.236138      TS             FALSE        TRUE
#> 9   2022-10-08 20:43:36     4.720493      TS             FALSE        TRUE
#> 15  2022-10-09 15:10:59     3.744255      TS             FALSE        TRUE
#> 17  2022-10-09 16:34:34     4.603939      TS             FALSE       FALSE
#> 24  2022-10-09 20:20:27     4.458669      TS             FALSE       FALSE
#> 32  2022-10-10 15:16:16     4.603857      TS             FALSE       FALSE
#> 146 2022-10-11 14:10:39     3.535534      TS             FALSE       FALSE
#> 152 2022-10-11 15:11:01     3.535534      TS             FALSE       FALSE
#> 165 2022-10-11 18:15:37     4.768179      TS             FALSE       FALSE
#>     tag:Early Intrusives, Intermediate-felsic (TTG)
#> 1                                             FALSE
#> 5                                             FALSE
#> 9                                              TRUE
#> 15                                            FALSE
#> 17                                            FALSE
#> 24                                            FALSE
#> 32                                            FALSE
#> 146                                           FALSE
#> 152                                           FALSE
#> 165                                           FALSE
#>     tag:Late Intrusives, felsic-mafic tag:Meta-Sedimentary rock, Quetico
#> 1                               FALSE                              FALSE
#> 5                               FALSE                              FALSE
#> 9                               FALSE                              FALSE
#> 15                              FALSE                              FALSE
#> 17                              FALSE                              FALSE
#> 24                              FALSE                              FALSE
#> 32                              FALSE                              FALSE
#> 146                             FALSE                              FALSE
#> 152                             FALSE                              FALSE
#> 165                             FALSE                              FALSE
#>     tag:Meta-Sedimentary rock, chemical tag:Metavolcanic, ultramafic-mafic
#> 1                                 FALSE                              FALSE
#> 5                                 FALSE                              FALSE
#> 9                                 FALSE                              FALSE
#> 15                                FALSE                              FALSE
#> 17                                FALSE                              FALSE
#> 24                                FALSE                              FALSE
#> 32                                FALSE                              FALSE
#> 146                               FALSE                              FALSE
#> 152                               FALSE                              FALSE
#> 165                               FALSE                              FALSE
#>     tag:Mylonitic schist
#> 1                  FALSE
#> 5                  FALSE
#> 9                  FALSE
#> 15                 FALSE
#> 17                 FALSE
#> 24                 FALSE
#> 32                 FALSE
#> 146                FALSE
#> 152                FALSE
#> 165                FALSE
#>     tag:SGSB Intrusives, Intermediate-felsic hypabyssal rocks
#> 1                                                        TRUE
#> 5                                                        TRUE
#> 9                                                        TRUE
#> 15                                                       TRUE
#> 17                                                      FALSE
#> 24                                                      FALSE
#> 32                                                      FALSE
#> 146                                                     FALSE
#> 152                                                     FALSE
#> 165                                                     FALSE
#>     tag:SGSB Intrusives, Mafic-Ultramafic
#> 1                                    TRUE
#> 5                                    TRUE
#> 9                                    TRUE
#> 15                                   TRUE
#> 17                                  FALSE
#> 24                                  FALSE
#> 32                                  FALSE
#> 146                                 FALSE
#> 152                                 FALSE
#> 165                                 FALSE
#>     tag:SGSB Volcanics, Intermediate-felsic
#> 1                                     FALSE
#> 5                                     FALSE
#> 9                                     FALSE
#> 15                                    FALSE
#> 17                                     TRUE
#> 24                                     TRUE
#> 32                                    FALSE
#> 146                                    TRUE
#> 152                                    TRUE
#> 165                                   FALSE
#>     tag:SGSB Volcanics, Intermediate-mafic                   geometry
#> 1                                    FALSE    POINT (-90.712 48.5426)
#> 5                                    FALSE POINT (-90.72518 48.53076)
#> 9                                    FALSE POINT (-90.72712 48.57296)
#> 15                                   FALSE POINT (-90.58906 48.60464)
#> 17                                   FALSE  POINT (-90.5811 48.60599)
#> 24                                   FALSE POINT (-90.57892 48.56819)
#> 32                                    TRUE POINT (-90.60366 48.54289)
#> 146                                  FALSE POINT (-90.74599 48.51911)
#> 152                                  FALSE POINT (-90.76746 48.52015)
#> 165                                   TRUE  POINT (-90.7947 48.53087)
#> 
#> $tags
#>          tag_type                                              tag_name
#>            <char>                                                <char>
#>  1: geologic_unit                                         Diabase dikes
#>  2: geologic_unit                                               Diorite
#>  3: geologic_unit           Early Intrusives, Intermediate-felsic (TTG)
#>  4: geologic_unit                         Late Intrusives, felsic-mafic
#>  5: geologic_unit                       Meta-Sedimentary rock, chemical
#>  6: geologic_unit                        Meta-Sedimentary rock, Quetico
#>  7: geologic_unit                           Meta-Sedimentary rock, Wawa
#>  8: geologic_unit                        Metavolcanic, ultramafic-mafic
#>  9: geologic_unit                                      Mylonitic schist
#> 10: geologic_unit SGSB Intrusives, Intermediate-felsic hypabyssal rocks
#> 11: geologic_unit                     SGSB Intrusives, Mafic-Ultramafic
#> 12: geologic_unit                   SGSB Volcanics, Intermediate-felsic
#> 13: geologic_unit                    SGSB Volcanics, Intermediate-mafic
#>     tag_unit_label_abbreviation tag_rock_type tag_igneous_rock_class
#>                          <char>        <char>                 <char>
#>  1:                   Intrusion       igneous               plutonic
#>  2:                     Diorite          <NA>                   <NA>
#>  3:                   Intrusion       igneous               plutonic
#>  4:              Late Intrusion       igneous               plutonic
#>  5:                        Wawa   sedimentary                   <NA>
#>  6:                     Quetico   sedimentary                   <NA>
#>  7:                        Wawa   sedimentary                   <NA>
#>  8:              Ultramafic MVR   metamorphic                   <NA>
#>  9:                         Myl          <NA>                   <NA>
#> 10:                   Intrusion       igneous             hypabyssal
#> 11:                   Intrusion       igneous               plutonic
#> 12:                 Shebandowan       igneous               volcanic
#> 13:                 Shebandowan       igneous               volcanic
#>     tag_plutonic_rock_types tag_other_plutonic_rock_type       tag_id tag_color
#>                      <char>                       <char>        <num>    <char>
#>  1:                   other                      Diabase 1.665001e+13   #FF8000
#>  2:                    <NA>                         <NA> 1.665088e+13   #4C0099
#>  3:                   other          Intermediate-felsic 1.664997e+13   #CC00CC
#>  4:                    <NA>                         <NA> 1.664997e+13   #CC0000
#>  5:                    <NA>                         <NA> 1.664999e+13   #663300
#>  6:                    <NA>                         <NA> 1.664998e+13   #C0C0C0
#>  7:                    <NA>                         <NA> 1.664997e+13   #606060
#>  8:                    <NA>                         <NA> 1.687833e+13   #0066CC
#>  9:                    <NA>                         <NA> 1.665780e+13   #000000
#> 10:                    <NA>                         <NA> 1.664998e+13   #FF3399
#> 11:                    <NA>                         <NA> 1.664998e+13   #00FFFF
#> 12:                    <NA>                         <NA> 1.664999e+13   #4C9900
#> 13:                    <NA>                         <NA> 1.664999e+13   #00CC66
#>                                    tag_notes
#>                                       <char>
#>  1:                                       10
#>  2:                                     <NA>
#>  3:                                        8
#>  4: 9; e.g. Hood Lake Stock, Moss Lake Stock
#>  5:                                        3
#>  6:                                        4
#>  7:                                        7
#>  8:                                     <NA>
#>  9:                                     <NA>
#> 10:                                        6
#> 11:                                        5
#> 12:                                        2
#> 13:                                        1
#>                                                                               tag_description
#>                                                                                        <char>
#>  1:                                                                                      <NA>
#>  2:                                                                                      <NA>
#>  3:                                      Tonalite, trondhjemite, granodiorite, Bt-Ms granites
#>  4:                                         Syenite, monzonite, granite, diorite, lamprophyre
#>  5:                                                                        Cherts, ironstones
#>  6:                   Wackes, siltstones, mudstones, phyllite, arenites, partially migmatitic
#>  7:                                                                             Conglomerates
#>  8:                                                          Massive flow, serpentinized rock
#>  9:                                                                                      <NA>
#> 10:                                  Fsp porphyry, Qtz-Fsp porphyry, Amph￼-Py-Qz-Fsp porphyry
#> 11: Gabbro, diorite, leucogabbro, Chl-Act-schist, Amphibolite, Pl-phyric gabbro, dikes, sills
#> 12:                                Tuffs, massive flows, Qz-chlorite schists, volcanoclastics
#> 13:                            Pillows, chlorite schists, tuffs, tuff breccias, amphibolites,
#>               tag_map_unit_name tag_sedimentary_rock_type
#>                          <char>                    <char>
#>  1:                        <NA>                      <NA>
#>  2:                        <NA>                      <NA>
#>  3:                        <NA>                      <NA>
#>  4:                        <NA>                      <NA>
#>  5:                        <NA>                      <NA>
#>  6:                     Quetico                      <NA>
#>  7:                Wawa-Abitibi              conglomerate
#>  8: Shebandowan greenstone belt                      <NA>
#>  9:                        <NA>                      <NA>
#> 10:                        <NA>                      <NA>
#> 11:                        <NA>                      <NA>
#> 12: Shebandowan greenstone belt                      <NA>
#> 13: Shebandowan greenstone belt                      <NA>
#>     tag_metamorphic_rock_types tag_absolute_age_of_geologic_unit
#>                         <char>                             <num>
#>  1:                       <NA>                                NA
#>  2:                       <NA>                                NA
#>  3:                       <NA>                                NA
#>  4:                       <NA>                                NA
#>  5:                       <NA>                                NA
#>  6:                       <NA>                                NA
#>  7:                       <NA>                                NA
#>  8:             meta_ultramafi                              2720
#>  9:                       <NA>                                NA
#> 10:                       <NA>                                NA
#> 11:                       <NA>                                NA
#> 12:                       <NA>                                NA
#> 13:                       <NA>                                NA
#>     tag_age_uncertainty tag_metamorphic_grade tag_volcanic_rock_type
#>                   <num>                <char>                 <char>
#>  1:                  NA                  <NA>                   <NA>
#>  2:                  NA                  <NA>                   <NA>
#>  3:                  NA                  <NA>                   <NA>
#>  4:                  NA                  <NA>                   <NA>
#>  5:                  NA                  <NA>                   <NA>
#>  6:                  NA                  <NA>                   <NA>
#>  7:                  NA                  <NA>                   <NA>
#>  8:                  10        greenschist_fa                   <NA>
#>  9:                  NA                  <NA>                   <NA>
#> 10:                  NA                  <NA>                   <NA>
#> 11:                  NA                  <NA>                   <NA>
#> 12:                  NA                  <NA>                   <NA>
#> 13:                  NA                  <NA>                  other
#>     tag_other_volcanic_rock_type
#>                           <char>
#>  1:                         <NA>
#>  2:                         <NA>
#>  3:                         <NA>
#>  4:                         <NA>
#>  5:                         <NA>
#>  6:                         <NA>
#>  7:                         <NA>
#>  8:                         <NA>
#>  9:                         <NA>
#> 10:                         <NA>
#> 11:                         <NA>
#> 12:                         <NA>
#> 13:           Intermediate-mafic
#> 
#> $planar
#> Plane object (n = 323):
#>        dip_direction dip
#>   [1,]           332  87
#>   [2,]           337  89
#>   [3,]           145  84
#>   [4,]           330  82
#>   [5,]           334  81
#>   [6,]           340  85
#>   [7,]           167  89
#>   [8,]           138  84
#>   [9,]           306  77
#>  [10,]           186  83
#>  [11,]           313  82
#>  [12,]           148  79
#>  [13,]           142  79
#>  [14,]           142  79
#>  [15,]           151  76
#>  [16,]           165  78
#>  [17,]           167  68
#>  [18,]           126  73
#>  [19,]           127  74
#>  [20,]           139  73
#>  [21,]           194  76
#>  [22,]           160  85
#>  [23,]           166  81
#>  [24,]           165  75
#>  [25,]           162  84
#>  [26,]           165  81
#>  [27,]           166  85
#>  [28,]           152  47
#>  [29,]           323  75
#>  [30,]           337  72
#>  [31,]           163  79
#>  [32,]           343  84
#>  [33,]           156  88
#>  [34,]           161  80
#>  [35,]           348  86
#>  [36,]           328  89
#>  [37,]           337  77
#>  [38,]           109  67
#>  [39,]           113  83
#>  [40,]           113  75
#>  [41,]           115  68
#>  [42,]           117  63
#>  [43,]           106  73
#>  [44,]           111  70
#>  [45,]           117  69
#>  [46,]           129  66
#>  [47,]           129  69
#>  [48,]           134  73
#>  [49,]           111  74
#>  [50,]           108  65
#>  [51,]           118  80
#>  [52,]           297  74
#>  [53,]           143  79
#>  [54,]           140  77
#>  [55,]           153  88
#>  [56,]           158  77
#>  [57,]           137  78
#>  [58,]           128  80
#>  [59,]           327  89
#>  [60,]           329  83
#>  [61,]           334  87
#>  [62,]           144  87
#>  [63,]           146  84
#>  [64,]           334  87
#>  [65,]           142  80
#>  [66,]           151  78
#>  [67,]           137  56
#>  [68,]           141  75
#>  [69,]           163  85
#>  [70,]           336  74
#>  [71,]           345  81
#>  [72,]           159  87
#>  [73,]           165  89
#>  [74,]           170  81
#>  [75,]           346  82
#>  [76,]           338  85
#>  [77,]           356  88
#>  [78,]           161  77
#>  [79,]           158  88
#>  [80,]           159  87
#>  [81,]           337  89
#>  [82,]           349  86
#>  [83,]           350  84
#>  [84,]           155  84
#>  [85,]           348  82
#>  [86,]           342  89
#>  [87,]           340  84
#>  [88,]           338  84
#>  [89,]           349  85
#>  [90,]           154  89
#>  [91,]           336  84
#>  [92,]           352  89
#>  [93,]           344  88
#>  [94,]           171  88
#>  [95,]           351  88
#>  [96,]           180  87
#>  [97,]           346  82
#>  [98,]           320  74
#>  [99,]           135  85
#> [100,]           309  85
#> [101,]           326  80
#> [102,]           331  70
#> [103,]           329  66
#> [104,]           322  74
#> [105,]           330  78
#> [106,]           311  78
#> [107,]           146  88
#> [108,]           332  87
#> [109,]           338  87
#> [110,]           186  78
#> [111,]             3  81
#> [112,]           195  82
#> [113,]           150  81
#> [114,]           311  72
#> [115,]           315  79
#> [116,]           333  85
#> [117,]           313  81
#> [118,]           162  76
#> [119,]           286  68
#> [120,]           161  87
#> [121,]           335  81
#> [122,]           327  83
#> [123,]           313  68
#> [124,]           316  71
#> [125,]           333  75
#> [126,]           319  89
#> [127,]           332  74
#> [128,]           149  89
#> [129,]            86  38
#> [130,]           320  68
#> [131,]           285  86
#> [132,]           330  76
#> [133,]           337  83
#> [134,]           130  89
#> [135,]           304  71
#> [136,]           148  86
#> [137,]           346  85
#> [138,]           347  89
#> [139,]           339  83
#> [140,]           151  84
#> [141,]           318  87
#> [142,]           146  89
#> [143,]           169  89
#> [144,]           328  88
#> [145,]           132  74
#> [146,]           321  47
#> [147,]           191  47
#> [148,]           299  48
#> [149,]           178  78
#> [150,]           191  84
#> [151,]           295  77
#> [152,]           305  89
#> [153,]           192  76
#> [154,]           164  73
#> [155,]           160  55
#> [156,]           358  83
#> [157,]           358  87
#> [158,]           336  85
#> [159,]           332  78
#> [160,]           164  85
#> [161,]           162  88
#> [162,]           331  87
#> [163,]           161  84
#> [164,]           169  77
#> [165,]           156  85
#> [166,]           146  79
#> [167,]           163  64
#> [168,]           128  85
#> [169,]           339  78
#> [170,]           338  77
#> [171,]           331  82
#> [172,]           332  84
#> [173,]           150  74
#> [174,]           341  87
#> [175,]           322  54
#> [176,]           299  70
#> [177,]           169  67
#> [178,]           166  75
#> [179,]           154  70
#> [180,]           353  89
#> [181,]           337  81
#> [182,]           330  87
#> [183,]           336  85
#> [184,]           193  79
#> [185,]           329  83
#> [186,]           330  84
#> [187,]           325  66
#> [188,]           324  81
#> [189,]           329  80
#> [190,]           344  64
#> [191,]           166  89
#> [192,]           339  85
#> [193,]           341  62
#> [194,]           338  58
#> [195,]           159  79
#> [196,]           341  80
#> [197,]           335  79
#> [198,]           345  50
#> [199,]           344  52
#> [200,]           346  60
#> [201,]           352  53
#> [202,]             6  42
#> [203,]             6  48
#> [204,]             4  38
#> [205,]             0  42
#> [206,]             2  41
#> [207,]           349  41
#> [208,]           352  36
#> [209,]           348  41
#> [210,]           348  52
#> [211,]           357  36
#> [212,]           346  65
#> [213,]           352  42
#> [214,]           354  30
#> [215,]             2  39
#> [216,]           346  46
#> [217,]           148  42
#> [218,]           150  46
#> [219,]           139  51
#> [220,]           158  51
#> [221,]           152  43
#> [222,]           150  48
#> [223,]           145  49
#> [224,]           146  63
#> [225,]           157  69
#> [226,]           343  74
#> [227,]           345  71
#> [228,]           349  74
#> [229,]           344  76
#> [230,]           345  74
#> [231,]           316  85
#> [232,]           145  86
#> [233,]           311  46
#> [234,]           297  40
#> [235,]           323  77
#> [236,]           318  73
#> [237,]           314  88
#> [238,]           319  71
#> [239,]           159  86
#> [240,]           335  84
#> [241,]           339  86
#> [242,]           330  89
#> [243,]           326  88
#> [244,]           334  81
#> [245,]           330  89
#> [246,]           333  81
#> [247,]           336  88
#> [248,]           332  84
#> [249,]           334  78
#> [250,]           141  88
#> [251,]           327  82
#> [252,]           341  81
#> [253,]           190  33
#> [254,]           195  36
#> [255,]           167  30
#> [256,]           194  32
#> [257,]           187  35
#> [258,]           342  86
#> [259,]           334  78
#> [260,]           353  82
#> [261,]             4  72
#> [262,]             5  32
#> [263,]             0  61
#> [264,]             9  73
#> [265,]            10  36
#> [266,]            23  35
#> [267,]           326  26
#> [268,]           328  37
#> [269,]           339  60
#> [270,]           338  58
#> [271,]             8  59
#> [272,]             5  58
#> [273,]           353  35
#> [274,]           329  66
#> [275,]           331  70
#> [276,]           337  74
#> [277,]           333  81
#> [278,]           157  82
#> [279,]            12  22
#> [280,]            89  57
#> [281,]            62  40
#> [282,]            34  64
#> [283,]           350  59
#> [284,]            14  51
#> [285,]           353  60
#> [286,]           119  56
#> [287,]           323  73
#> [288,]           341  89
#> [289,]           350  29
#> [290,]           348  50
#> [291,]           334  55
#> [292,]           342  86
#> [293,]           350  74
#> [294,]           345  71
#> [295,]           329  73
#> [296,]           336  72
#> [297,]           341  73
#> [298,]           306  59
#> [299,]           331  66
#> [300,]           168  68
#> [301,]           161  72
#> [302,]           348  86
#> [303,]           162  85
#> [304,]           173  89
#> [305,]           342  78
#> [306,]           164  87
#> [307,]           343  74
#> [308,]           161  88
#> [309,]            50  60
#> [310,]           342  57
#> [311,]           347  82
#> [312,]           210  77
#> [313,]           197  59
#> [314,]           142  65
#> [315,]           141  68
#> [316,]           332  88
#> [317,]           140  83
#> [318,]           348  88
#> [319,]           325  82
#> [320,]           314  73
#> [321,]           158  90
#> [322,]           347  88
#> [323,]           320  82
#> 
#> $linear
#> Line object (n = 323):
#>        azimuth plunge
#>   [1,]      62      8
#>   [2,]      67      4
#>   [3,]     235      2
#>   [4,]      60      0
#>   [5,]     246     14
#>   [6,]      69     13
#>   [7,]      77      5
#>   [8,]      48      4
#>   [9,]      34      7
#>  [10,]      96      4
#>  [11,]      37     38
#>  [12,]     228     44
#>  [13,]     220     48
#>  [14,]     224     37
#>  [15,]     208     66
#>  [16,]     199     76
#>  [17,]      87     24
#>  [18,]      49     36
#>  [19,]      53     46
#>  [20,]      72     53
#>  [21,]     105      4
#>  [22,]      73     27
#>  [23,]      76      4
#>  [24,]      81     20
#>  [25,]      72      1
#>  [26,]     255      1
#>  [27,]      77     10
#>  [28,]     241      1
#>  [29,]      48     17
#>  [30,]      62     16
#>  [31,]      78     25
#>  [32,]      70     22
#>  [33,]      67     22
#>  [34,]      74     20
#>  [35,]      76     23
#>  [36,]      58      4
#>  [37,]      62     24
#>  [38,]     183     35
#>  [39,]      24      7
#>  [40,]      26     11
#>  [41,]      29     10
#>  [42,]      31      7
#>  [43,]      17      1
#>  [44,]      22      3
#>  [45,]      28      2
#>  [46,]     213     13
#>  [47,]     213     15
#>  [48,]     217     23
#>  [49,]     201      0
#>  [50,]      19      3
#>  [51,]      29      7
#>  [52,]      20     24
#>  [53,]      58     24
#>  [54,]      55     22
#>  [55,]      63     12
#>  [56,]     247      2
#>  [57,]      47      3
#>  [58,]      39      3
#>  [59,]      56      4
#>  [60,]     239      0
#>  [61,]      64     10
#>  [62,]      55     12
#>  [63,]      57      7
#>  [64,]      64      5
#>  [65,]      55     16
#>  [66,]      62      4
#>  [67,]      52      7
#>  [68,]      52      4
#>  [69,]      75     19
#>  [70,]     285     66
#>  [71,]      73      9
#>  [72,]      70      6
#>  [73,]     255     13
#>  [74,]      82     15
#>  [75,]      70     38
#>  [76,]      61     55
#>  [77,]     266      0
#>  [78,]      74     15
#>  [79,]      69     38
#>  [80,]      70     22
#>  [81,]      67     16
#>  [82,]     259      0
#>  [83,]      76     31
#>  [84,]      67     21
#>  [85,]      75     17
#>  [86,]      71     12
#>  [87,]      64     46
#>  [88,]      65     23
#>  [89,]      76     32
#>  [90,]     244      4
#>  [91,]      63     23
#>  [92,]      82     17
#>  [93,]      74     23
#>  [94,]     261      0
#>  [95,]     261      6
#>  [96,]      91     27
#>  [97,]      75     10
#>  [98,]     231      6
#>  [99,]      48     33
#> [100,]      38      7
#> [101,]      51     30
#> [102,]      58     10
#> [103,]      58      1
#> [104,]      48     15
#> [105,]      57     15
#> [106,]      29     43
#> [107,]      57     19
#> [108,]      62      8
#> [109,]     248     13
#> [110,]      98      6
#> [111,]      92      4
#> [112,]     106      2
#> [113,]      66     33
#> [114,]      41      0
#> [115,]     227     10
#> [116,]     243      8
#> [117,]     223      0
#> [118,]     252      0
#> [119,]     198      5
#> [120,]      71     10
#> [121,]      62     18
#> [122,]     237      1
#> [123,]      35     17
#> [124,]     226      0
#> [125,]      60     11
#> [126,]      49     10
#> [127,]      56     21
#> [128,]     239     14
#> [129,]      94     37
#> [130,]     236     13
#> [131,]      14     17
#> [132,]     245     21
#> [133,]     250     25
#> [134,]      40     19
#> [135,]      31      8
#> [136,]      60     35
#> [137,]      75     15
#> [138,]      77     16
#> [139,]      66     23
#> [140,]      64     29
#> [141,]     228      1
#> [142,]     236     19
#> [143,]      79      0
#> [144,]      57     21
#> [145,]     163     71
#> [146,]     268     33
#> [147,]     227     41
#> [148,]     252     37
#> [149,]     266     10
#> [150,]     103     19
#> [151,]     206      4
#> [152,]      35     73
#> [153,]     281      4
#> [154,]      92     46
#> [155,]      95     33
#> [156,]      87     10
#> [157,]     269     14
#> [158,]      65      8
#> [159,]      60     13
#> [160,]      76     24
#> [161,]      72     16
#> [162,]      60      9
#> [163,]      74     26
#> [164,]      83     16
#> [165,]      67      9
#> [166,]      57      5
#> [167,]      81     15
#> [168,]      39     10
#> [169,]      67      8
#> [170,]      67      4
#> [171,]      59     16
#> [172,]      61     11
#> [173,]      66     21
#> [174,]      71      3
#> [175,]      16     38
#> [176,]     213     12
#> [177,]      86     18
#> [178,]      84     28
#> [179,]      79     35
#> [180,]      83     17
#> [181,]      66      1
#> [182,]     240     10
#> [183,]     248     22
#> [184,]     104      7
#> [185,]     240     10
#> [186,]     241     12
#> [187,]     236      2
#> [188,]     234      0
#> [189,]     240      2
#> [190,]      61     25
#> [191,]      76      4
#> [192,]      67     17
#> [193,]      62     16
#> [194,]      61     11
#> [195,]      71      8
#> [196,]      70      8
#> [197,]      63     13
#> [198,]      63     13
#> [199,]      67      9
#> [200,]      71      9
#> [201,]      73     11
#> [202,]     291     13
#> [203,]      79     17
#> [204,]     278      3
#> [205,]      83      6
#> [206,]      79     11
#> [207,]      72      5
#> [208,]      64     12
#> [209,]      68      8
#> [210,]      61     20
#> [211,]      71     11
#> [212,]      68     16
#> [213,]      68     12
#> [214,]      69      8
#> [215,]      71     16
#> [216,]      52     22
#> [217,]      74     14
#> [218,]      80     19
#> [219,]      63     17
#> [220,]      76     10
#> [221,]      81     16
#> [222,]      75     16
#> [223,]      67     13
#> [224,]      64     15
#> [225,]      70      7
#> [226,]     254      3
#> [227,]     258      9
#> [228,]      77      5
#> [229,]     254      0
#> [230,]     256      1
#> [231,]     226      0
#> [232,]      55      9
#> [233,]     228      7
#> [234,]     222     12
#> [235,]      50     15
#> [236,]      47      4
#> [237,]      43     22
#> [238,]      41     23
#> [239,]      70     11
#> [240,]     248     30
#> [241,]     253     39
#> [242,]     240     14
#> [243,]     237     16
#> [244,]     252     46
#> [245,]     240     35
#> [246,]     249     30
#> [247,]     247     38
#> [248,]     245     31
#> [249,]     245      3
#> [250,]     230     44
#> [251,]      57      0
#> [252,]     251      5
#> [253,]     261     12
#> [254,]     197     36
#> [255,]     174     30
#> [256,]     200     31
#> [257,]     151     30
#> [258,]      72      1
#> [259,]      39     64
#> [260,]      72     54
#> [261,]      65     57
#> [262,]     347     31
#> [263,]     312     50
#> [264,]      93     17
#> [265,]      61     25
#> [266,]      65     28
#> [267,]      54      1
#> [268,]      51      5
#> [269,]     252      5
#> [270,]     255     10
#> [271,]      86     20
#> [272,]      81     22
#> [273,]      64     13
#> [274,]     245     11
#> [275,]     244      8
#> [276,]      65      6
#> [277,]      63      3
#> [278,]      67      1
#> [279,]      31     20
#> [280,]      73     56
#> [281,]      49     40
#> [282,]     111     25
#> [283,]      66     22
#> [284,]     293     10
#> [285,]      71     20
#> [286,]      51     29
#> [287,]      36     45
#> [288,]      70     25
#> [289,]      67      7
#> [290,]      60     20
#> [291,]      53     16
#> [292,]      71     12
#> [293,]      74     21
#> [294,]     261     17
#> [295,]      49     30
#> [296,]      61     13
#> [297,]      57     39
#> [298,]     218      2
#> [299,]      54     15
#> [300,]     111     53
#> [301,]      76     14
#> [302,]      78      0
#> [303,]      77     48
#> [304,]      83      8
#> [305,]     253      7
#> [306,]     254     12
#> [307,]      72      3
#> [308,]      71      4
#> [309,]     139      2
#> [310,]      61     18
#> [311,]      76      4
#> [312,]     125     22
#> [313,]     124     26
#> [314,]     194     53
#> [315,]     225     13
#> [316,]      62      8
#> [317,]     230      2
#> [318,]      78      9
#> [319,]     237     12
#> [320,]      44      1
#> [321,]     248     15
#> [322,]     258      4
#> [323,]     232     14
#> 
```
