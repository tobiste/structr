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
#>                                        id dip_direction   dip strike plunge
#>                                    <char>         <num> <num>  <num>  <num>
#>   1: d951938c-6bee-4ce4-bee7-6cf17693a090           332    87    242      8
#>   2: 4a0aba4f-4707-4852-9811-f8457fd226e7           337    89    247      4
#>   3: 86c34f8d-b277-4a27-8f49-9339df976f66           145    84     55      2
#>   4: ec15d95a-e23c-4ba5-8280-63019e7ce553           330    82    240      0
#>   5: ad910789-d6b2-4179-9c38-9e95149d7b56           334    81    244     14
#>  ---                                                                       
#> 343: 3ae2760c-136b-4cb4-bc37-95c852642185           325    82    235     12
#> 344: 2810a725-9eee-47d7-9c97-5b5812092240           314    73    224      1
#> 345: 88d22c99-21ad-460d-b364-a8bd6dab7c7e           158    90     68     15
#> 346: 64086f81-b490-40f7-8291-f37d8c903d41           347    88    257      4
#> 347: e9b4bca6-37d7-4a2b-8833-ddc711da0a79           320    82    230     14
#>      trend associated planar_type quality unix_timestamp  notes
#>      <num>     <lgcl>      <char>  <char>          <num> <char>
#>   1:    62       TRUE        <NA>       5   1.665089e+12   <NA>
#>   2:    67       TRUE        <NA>       5   1.665089e+12   <NA>
#>   3:   235       TRUE        <NA>       5   1.665089e+12   <NA>
#>   4:    60       TRUE        <NA>       5   1.665089e+12   <NA>
#>   5:   246       TRUE        <NA>       5   1.665089e+12   <NA>
#>  ---                                                           
#> 343:   237       TRUE        <NA>       4   1.721928e+12   <NA>
#> 344:    44       TRUE        <NA>       5   1.721928e+12   <NA>
#> 345:   248       TRUE        <NA>       3   1.721928e+12   <NA>
#> 346:   258       TRUE        <NA>       3   1.721928e+12   <NA>
#> 347:   232       TRUE        <NA>       5   1.721930e+12   <NA>
#>      modified_timestamp          spot        spot_id feature_type
#>                   <num>        <char>         <char>       <char>
#>   1:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   2:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   3:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   4:                 NA 22-TS-Moss-02 16650869118218    foliation
#>   5:                 NA 22-TS-Moss-02 16650869118218    foliation
#>  ---                                                             
#> 343:                 NA       24TS-16 17219277194544    foliation
#> 344:                 NA       24TS-16 17219277194544    foliation
#> 345:                 NA       24TS-16 17219277194544    foliation
#> 346:                 NA       24TS-16 17219277194544    foliation
#> 347:                 NA       24TS-17 17219297078328    foliation
#>      foliation_type         label contact_type fault_or_sz_type bedding_type
#>              <char>        <char>       <char>           <char>       <char>
#>   1:       cleavage foliation 242         <NA>             <NA>         <NA>
#>   2:       cleavage foliation 247         <NA>             <NA>         <NA>
#>   3:       cleavage  foliation 55         <NA>             <NA>         <NA>
#>   4:       cleavage foliation 240         <NA>             <NA>         <NA>
#>   5:       cleavage foliation 244         <NA>             <NA>         <NA>
#>  ---                                                                        
#> 343:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 344:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 345:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 346:           <NA>          <NA>         <NA>             <NA>         <NA>
#> 347:           <NA>          <NA>         <NA>             <NA>         <NA>
#>      vein_fill directional_indicators vein_type other_vein_fill
#>         <char>                 <char>    <char>          <char>
#>   1:      <NA>                   <NA>      <NA>            <NA>
#>   2:      <NA>                   <NA>      <NA>            <NA>
#>   3:      <NA>                   <NA>      <NA>            <NA>
#>   4:      <NA>                   <NA>      <NA>            <NA>
#>   5:      <NA>                   <NA>      <NA>            <NA>
#>  ---                                                           
#> 343:      <NA>                   <NA>      <NA>            <NA>
#> 344:      <NA>                   <NA>      <NA>            <NA>
#> 345:      <NA>                   <NA>      <NA>            <NA>
#> 346:      <NA>                   <NA>      <NA>            <NA>
#> 347:      <NA>                   <NA>      <NA>            <NA>
#>      foliation_defined_by movement other_feature fracture_type
#>                    <char>   <char>        <char>        <char>
#>   1:                 <NA>     <NA>          <NA>          <NA>
#>   2:                 <NA>     <NA>          <NA>          <NA>
#>   3:                 <NA>     <NA>          <NA>          <NA>
#>   4:                 <NA>     <NA>          <NA>          <NA>
#>   5:                 <NA>     <NA>          <NA>          <NA>
#>  ---                                                          
#> 343:                 <NA>     <NA>          <NA>          <NA>
#> 344:                 <NA>     <NA>          <NA>          <NA>
#> 345:                 <NA>     <NA>          <NA>          <NA>
#> 346:                 <NA>     <NA>          <NA>          <NA>
#> 347:                 <NA>     <NA>          <NA>          <NA>
#>      movement_justification        linear_type linear_quality linear_notes
#>                      <char>             <char>         <char>       <char>
#>   1:                   <NA> linear_orientation              5         <NA>
#>   2:                   <NA> linear_orientation              5         <NA>
#>   3:                   <NA> linear_orientation              5         <NA>
#>   4:                   <NA> linear_orientation              5         <NA>
#>   5:                   <NA> linear_orientation              5         <NA>
#>  ---                                                                      
#> 343:                   <NA> linear_orientation              4         <NA>
#> 344:                   <NA> linear_orientation              5         <NA>
#> 345:                   <NA> linear_orientation              3         <NA>
#> 346:                   <NA> linear_orientation              3         <NA>
#> 347:                   <NA> linear_orientation              5         <NA>
#>      linear_feature_type linear_defined_by   linear_label               type
#>                   <char>            <char>         <char>             <char>
#>   1:          stretching              <NA>  stretching 62 planar_orientation
#>   2:          stretching              <NA>  stretching 67 planar_orientation
#>   3:          stretching              <NA> stretching 235 planar_orientation
#>   4:          stretching              <NA>  stretching 60 planar_orientation
#>   5:          stretching              <NA> stretching 246 planar_orientation
#>  ---                                                                        
#> 343:       mineral_align              <NA>           <NA> planar_orientation
#> 344:       mineral_align              <NA>           <NA> planar_orientation
#> 345:       mineral_align              <NA>           <NA> planar_orientation
#> 346:       mineral_align              <NA>           <NA> planar_orientation
#> 347:       mineral_align              <NA>           <NA> planar_orientation
#>      linear_unix_timestamp                            linear_id defined_by
#>                      <num>                               <char>     <char>
#>   1:          1.665089e+12 50b1dbef-5d64-4a8c-b427-7b1bb8c6053e       <NA>
#>   2:          1.665089e+12 af52e337-9c6b-4392-9f75-eb2dea09e016       <NA>
#>   3:          1.665089e+12 8189eda4-3acd-4bb8-9075-7cf5b6e3b2c1       <NA>
#>   4:          1.665089e+12 c8d4b893-6d70-4f89-b46b-ace81a01136d       <NA>
#>   5:          1.665089e+12 bf08e531-8c0e-4337-b9a3-7eccd8019486       <NA>
#>  ---                                                                      
#> 343:          1.721928e+12 84bcb839-a8fa-4a25-b28e-16b97ecbdce6       <NA>
#> 344:          1.721928e+12 9470b29b-550a-4613-9834-8667684f2aa6       <NA>
#> 345:          1.721928e+12 ceb03667-90cd-45bd-adfe-24102cbe7e42       <NA>
#> 346:          1.721928e+12 e55af915-d9a1-4d38-a710-71172ef34455       <NA>
#> 347:          1.721930e+12 9ff0f807-d54e-4dd9-b79c-57f1a1ff1116       <NA>
#>      linear_other_feature
#>                    <char>
#>   1:                 <NA>
#>   2:                 <NA>
#>   3:                 <NA>
#>   4:                 <NA>
#>   5:                 <NA>
#>  ---                     
#> 343:                 <NA>
#> 344:                 <NA>
#> 345:                 <NA>
#> 346:                 <NA>
#> 347:                 <NA>
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
#> Plane object (n = 347):
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
#>  [10,]           306  77
#>  [11,]           186  83
#>  [12,]           313  82
#>  [13,]           148  79
#>  [14,]           142  79
#>  [15,]           142  79
#>  [16,]           151  76
#>  [17,]           165  78
#>  [18,]           167  68
#>  [19,]           126  73
#>  [20,]           127  74
#>  [21,]           139  73
#>  [22,]           194  76
#>  [23,]           160  85
#>  [24,]           166  81
#>  [25,]           165  75
#>  [26,]           162  84
#>  [27,]           165  81
#>  [28,]           166  85
#>  [29,]           152  47
#>  [30,]           323  75
#>  [31,]           337  72
#>  [32,]           163  79
#>  [33,]           343  84
#>  [34,]           156  88
#>  [35,]           161  80
#>  [36,]           348  86
#>  [37,]           328  89
#>  [38,]           337  77
#>  [39,]           109  67
#>  [40,]           113  83
#>  [41,]           113  83
#>  [42,]           113  75
#>  [43,]           113  75
#>  [44,]           115  68
#>  [45,]           115  68
#>  [46,]           117  63
#>  [47,]           117  63
#>  [48,]           106  73
#>  [49,]           111  70
#>  [50,]           117  69
#>  [51,]           129  66
#>  [52,]           129  66
#>  [53,]           129  69
#>  [54,]           129  69
#>  [55,]           134  73
#>  [56,]           134  73
#>  [57,]           111  74
#>  [58,]           111  74
#>  [59,]           108  65
#>  [60,]           108  65
#>  [61,]           118  80
#>  [62,]           297  74
#>  [63,]           143  79
#>  [64,]           140  77
#>  [65,]           153  88
#>  [66,]           153  88
#>  [67,]           158  77
#>  [68,]           158  77
#>  [69,]           137  78
#>  [70,]           137  78
#>  [71,]           128  80
#>  [72,]           128  80
#>  [73,]           327  89
#>  [74,]           327  89
#>  [75,]           329  83
#>  [76,]           329  83
#>  [77,]           334  87
#>  [78,]           334  87
#>  [79,]           144  87
#>  [80,]           144  87
#>  [81,]           146  84
#>  [82,]           146  84
#>  [83,]           334  87
#>  [84,]           334  87
#>  [85,]           142  80
#>  [86,]           142  80
#>  [87,]           151  78
#>  [88,]           151  78
#>  [89,]           137  56
#>  [90,]           137  56
#>  [91,]           141  75
#>  [92,]           141  75
#>  [93,]           163  85
#>  [94,]           336  74
#>  [95,]           345  81
#>  [96,]           159  87
#>  [97,]           165  89
#>  [98,]           170  81
#>  [99,]           346  82
#> [100,]           338  85
#> [101,]           356  88
#> [102,]           161  77
#> [103,]           158  88
#> [104,]           159  87
#> [105,]           337  89
#> [106,]           349  86
#> [107,]           350  84
#> [108,]           155  84
#> [109,]           348  82
#> [110,]           342  89
#> [111,]           340  84
#> [112,]           338  84
#> [113,]           349  85
#> [114,]           154  89
#> [115,]           336  84
#> [116,]           352  89
#> [117,]           344  88
#> [118,]           171  88
#> [119,]           351  88
#> [120,]           180  87
#> [121,]           346  82
#> [122,]           320  74
#> [123,]           135  85
#> [124,]           309  85
#> [125,]           326  80
#> [126,]           331  70
#> [127,]           329  66
#> [128,]           322  74
#> [129,]           330  78
#> [130,]           311  78
#> [131,]           146  88
#> [132,]           332  87
#> [133,]           338  87
#> [134,]           186  78
#> [135,]             3  81
#> [136,]           195  82
#> [137,]           150  81
#> [138,]           311  72
#> [139,]           315  79
#> [140,]           333  85
#> [141,]           313  81
#> [142,]           162  76
#> [143,]           286  68
#> [144,]           161  87
#> [145,]           335  81
#> [146,]           327  83
#> [147,]           313  68
#> [148,]           316  71
#> [149,]           333  75
#> [150,]           319  89
#> [151,]           332  74
#> [152,]           149  89
#> [153,]            86  38
#> [154,]           320  68
#> [155,]           285  86
#> [156,]           330  76
#> [157,]           337  83
#> [158,]           130  89
#> [159,]           304  71
#> [160,]           148  86
#> [161,]           346  85
#> [162,]           347  89
#> [163,]           339  83
#> [164,]           151  84
#> [165,]           318  87
#> [166,]           146  89
#> [167,]           169  89
#> [168,]           328  88
#> [169,]           132  74
#> [170,]           321  47
#> [171,]           191  47
#> [172,]           299  48
#> [173,]           178  78
#> [174,]           191  84
#> [175,]           295  77
#> [176,]           305  89
#> [177,]           192  76
#> [178,]           164  73
#> [179,]           160  55
#> [180,]           358  83
#> [181,]           358  87
#> [182,]           336  85
#> [183,]           332  78
#> [184,]           164  85
#> [185,]           162  88
#> [186,]           331  87
#> [187,]           161  84
#> [188,]           169  77
#> [189,]           156  85
#> [190,]           146  79
#> [191,]           163  64
#> [192,]           128  85
#> [193,]           339  78
#> [194,]           338  77
#> [195,]           331  82
#> [196,]           332  84
#> [197,]           150  74
#> [198,]           341  87
#> [199,]           322  54
#> [200,]           299  70
#> [201,]           169  67
#> [202,]           166  75
#> [203,]           154  70
#> [204,]           353  89
#> [205,]           337  81
#> [206,]           330  87
#> [207,]           336  85
#> [208,]           193  79
#> [209,]           329  83
#> [210,]           330  84
#> [211,]           325  66
#> [212,]           324  81
#> [213,]           329  80
#> [214,]           344  64
#> [215,]           166  89
#> [216,]           339  85
#> [217,]           341  62
#> [218,]           338  58
#> [219,]           159  79
#> [220,]           341  80
#> [221,]           335  79
#> [222,]           345  50
#> [223,]           344  52
#> [224,]           346  60
#> [225,]           352  53
#> [226,]             6  42
#> [227,]             6  48
#> [228,]             4  38
#> [229,]             0  42
#> [230,]             2  41
#> [231,]           349  41
#> [232,]           352  36
#> [233,]           348  41
#> [234,]           348  52
#> [235,]           357  36
#> [236,]           346  65
#> [237,]           352  42
#> [238,]           354  30
#> [239,]             2  39
#> [240,]           346  46
#> [241,]           148  42
#> [242,]           150  46
#> [243,]           139  51
#> [244,]           158  51
#> [245,]           152  43
#> [246,]           150  48
#> [247,]           145  49
#> [248,]           146  63
#> [249,]           157  69
#> [250,]           343  74
#> [251,]           345  71
#> [252,]           349  74
#> [253,]           344  76
#> [254,]           345  74
#> [255,]           316  85
#> [256,]           145  86
#> [257,]           311  46
#> [258,]           297  40
#> [259,]           323  77
#> [260,]           318  73
#> [261,]           314  88
#> [262,]           319  71
#> [263,]           159  86
#> [264,]           335  84
#> [265,]           339  86
#> [266,]           330  89
#> [267,]           326  88
#> [268,]           334  81
#> [269,]           330  89
#> [270,]           333  81
#> [271,]           336  88
#> [272,]           332  84
#> [273,]           334  78
#> [274,]           141  88
#> [275,]           327  82
#> [276,]           341  81
#> [277,]           190  33
#> [278,]           195  36
#> [279,]           167  30
#> [280,]           194  32
#> [281,]           187  35
#> [282,]           342  86
#> [283,]           334  78
#> [284,]           353  82
#> [285,]             4  72
#> [286,]             5  32
#> [287,]             0  61
#> [288,]             9  73
#> [289,]            10  36
#> [290,]            23  35
#> [291,]           326  26
#> [292,]           328  37
#> [293,]           339  60
#> [294,]           338  58
#> [295,]             8  59
#> [296,]             5  58
#> [297,]           353  35
#> [298,]           329  66
#> [299,]           331  70
#> [300,]           337  74
#> [301,]           333  81
#> [302,]           157  82
#> [303,]            12  22
#> [304,]            89  57
#> [305,]            62  40
#> [306,]            34  64
#> [307,]           350  59
#> [308,]            14  51
#> [309,]           353  60
#> [310,]           119  56
#> [311,]           323  73
#> [312,]           341  89
#> [313,]           350  29
#> [314,]           348  50
#> [315,]           334  55
#> [316,]           342  86
#> [317,]           350  74
#> [318,]           345  71
#> [319,]           329  73
#> [320,]           336  72
#> [321,]           341  73
#> [322,]           306  59
#> [323,]           331  66
#> [324,]           168  68
#> [325,]           161  72
#> [326,]           348  86
#> [327,]           162  85
#> [328,]           173  89
#> [329,]           342  78
#> [330,]           164  87
#> [331,]           343  74
#> [332,]           161  88
#> [333,]            50  60
#> [334,]           342  57
#> [335,]           347  82
#> [336,]           210  77
#> [337,]           197  59
#> [338,]           142  65
#> [339,]           141  68
#> [340,]           332  88
#> [341,]           140  83
#> [342,]           348  88
#> [343,]           325  82
#> [344,]           314  73
#> [345,]           158  90
#> [346,]           347  88
#> [347,]           320  82
#> 
#> $linear
#> Line object (n = 347):
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
#>  [10,]      34      7
#>  [11,]      96      4
#>  [12,]      37     38
#>  [13,]     228     44
#>  [14,]     220     48
#>  [15,]     224     37
#>  [16,]     208     66
#>  [17,]     199     76
#>  [18,]      87     24
#>  [19,]      49     36
#>  [20,]      53     46
#>  [21,]      72     53
#>  [22,]     105      4
#>  [23,]      73     27
#>  [24,]      76      4
#>  [25,]      81     20
#>  [26,]      72      1
#>  [27,]     255      1
#>  [28,]      77     10
#>  [29,]     241      1
#>  [30,]      48     17
#>  [31,]      62     16
#>  [32,]      78     25
#>  [33,]      70     22
#>  [34,]      67     22
#>  [35,]      74     20
#>  [36,]      76     23
#>  [37,]      58      4
#>  [38,]      62     24
#>  [39,]     183     35
#>  [40,]      24      7
#>  [41,]      24      7
#>  [42,]      26     11
#>  [43,]      26     11
#>  [44,]      29     10
#>  [45,]      29     10
#>  [46,]      31      7
#>  [47,]      31      7
#>  [48,]      17      1
#>  [49,]      22      3
#>  [50,]      28      2
#>  [51,]     213     13
#>  [52,]     213     13
#>  [53,]     213     15
#>  [54,]     213     15
#>  [55,]     217     23
#>  [56,]     217     23
#>  [57,]     201      0
#>  [58,]     201      0
#>  [59,]      19      3
#>  [60,]      19      3
#>  [61,]      29      7
#>  [62,]      20     24
#>  [63,]      58     24
#>  [64,]      55     22
#>  [65,]      63     12
#>  [66,]      63     12
#>  [67,]     247      2
#>  [68,]     247      2
#>  [69,]      47      3
#>  [70,]      47      3
#>  [71,]      39      3
#>  [72,]      39      3
#>  [73,]      56      4
#>  [74,]      56      4
#>  [75,]     239      0
#>  [76,]     239      0
#>  [77,]      64     10
#>  [78,]      64     10
#>  [79,]      55     12
#>  [80,]      55     12
#>  [81,]      57      7
#>  [82,]      57      7
#>  [83,]      64      5
#>  [84,]      64      5
#>  [85,]      55     16
#>  [86,]      55     16
#>  [87,]      62      4
#>  [88,]      62      4
#>  [89,]      52      7
#>  [90,]      52      7
#>  [91,]      52      4
#>  [92,]      52      4
#>  [93,]      75     19
#>  [94,]     285     66
#>  [95,]      73      9
#>  [96,]      70      6
#>  [97,]     255     13
#>  [98,]      82     15
#>  [99,]      70     38
#> [100,]      61     55
#> [101,]     266      0
#> [102,]      74     15
#> [103,]      69     38
#> [104,]      70     22
#> [105,]      67     16
#> [106,]     259      0
#> [107,]      76     31
#> [108,]      67     21
#> [109,]      75     17
#> [110,]      71     12
#> [111,]      64     46
#> [112,]      65     23
#> [113,]      76     32
#> [114,]     244      4
#> [115,]      63     23
#> [116,]      82     17
#> [117,]      74     23
#> [118,]     261      0
#> [119,]     261      6
#> [120,]      91     27
#> [121,]      75     10
#> [122,]     231      6
#> [123,]      48     33
#> [124,]      38      7
#> [125,]      51     30
#> [126,]      58     10
#> [127,]      58      1
#> [128,]      48     15
#> [129,]      57     15
#> [130,]      29     43
#> [131,]      57     19
#> [132,]      62      8
#> [133,]     248     13
#> [134,]      98      6
#> [135,]      92      4
#> [136,]     106      2
#> [137,]      66     33
#> [138,]      41      0
#> [139,]     227     10
#> [140,]     243      8
#> [141,]     223      0
#> [142,]     252      0
#> [143,]     198      5
#> [144,]      71     10
#> [145,]      62     18
#> [146,]     237      1
#> [147,]      35     17
#> [148,]     226      0
#> [149,]      60     11
#> [150,]      49     10
#> [151,]      56     21
#> [152,]     239     14
#> [153,]      94     37
#> [154,]     236     13
#> [155,]      14     17
#> [156,]     245     21
#> [157,]     250     25
#> [158,]      40     19
#> [159,]      31      8
#> [160,]      60     35
#> [161,]      75     15
#> [162,]      77     16
#> [163,]      66     23
#> [164,]      64     29
#> [165,]     228      1
#> [166,]     236     19
#> [167,]      79      0
#> [168,]      57     21
#> [169,]     163     71
#> [170,]     268     33
#> [171,]     227     41
#> [172,]     252     37
#> [173,]     266     10
#> [174,]     103     19
#> [175,]     206      4
#> [176,]      35     73
#> [177,]     281      4
#> [178,]      92     46
#> [179,]      95     33
#> [180,]      87     10
#> [181,]     269     14
#> [182,]      65      8
#> [183,]      60     13
#> [184,]      76     24
#> [185,]      72     16
#> [186,]      60      9
#> [187,]      74     26
#> [188,]      83     16
#> [189,]      67      9
#> [190,]      57      5
#> [191,]      81     15
#> [192,]      39     10
#> [193,]      67      8
#> [194,]      67      4
#> [195,]      59     16
#> [196,]      61     11
#> [197,]      66     21
#> [198,]      71      3
#> [199,]      16     38
#> [200,]     213     12
#> [201,]      86     18
#> [202,]      84     28
#> [203,]      79     35
#> [204,]      83     17
#> [205,]      66      1
#> [206,]     240     10
#> [207,]     248     22
#> [208,]     104      7
#> [209,]     240     10
#> [210,]     241     12
#> [211,]     236      2
#> [212,]     234      0
#> [213,]     240      2
#> [214,]      61     25
#> [215,]      76      4
#> [216,]      67     17
#> [217,]      62     16
#> [218,]      61     11
#> [219,]      71      8
#> [220,]      70      8
#> [221,]      63     13
#> [222,]      63     13
#> [223,]      67      9
#> [224,]      71      9
#> [225,]      73     11
#> [226,]     291     13
#> [227,]      79     17
#> [228,]     278      3
#> [229,]      83      6
#> [230,]      79     11
#> [231,]      72      5
#> [232,]      64     12
#> [233,]      68      8
#> [234,]      61     20
#> [235,]      71     11
#> [236,]      68     16
#> [237,]      68     12
#> [238,]      69      8
#> [239,]      71     16
#> [240,]      52     22
#> [241,]      74     14
#> [242,]      80     19
#> [243,]      63     17
#> [244,]      76     10
#> [245,]      81     16
#> [246,]      75     16
#> [247,]      67     13
#> [248,]      64     15
#> [249,]      70      7
#> [250,]     254      3
#> [251,]     258      9
#> [252,]      77      5
#> [253,]     254      0
#> [254,]     256      1
#> [255,]     226      0
#> [256,]      55      9
#> [257,]     228      7
#> [258,]     222     12
#> [259,]      50     15
#> [260,]      47      4
#> [261,]      43     22
#> [262,]      41     23
#> [263,]      70     11
#> [264,]     248     30
#> [265,]     253     39
#> [266,]     240     14
#> [267,]     237     16
#> [268,]     252     46
#> [269,]     240     35
#> [270,]     249     30
#> [271,]     247     38
#> [272,]     245     31
#> [273,]     245      3
#> [274,]     230     44
#> [275,]      57      0
#> [276,]     251      5
#> [277,]     261     12
#> [278,]     197     36
#> [279,]     174     30
#> [280,]     200     31
#> [281,]     151     30
#> [282,]      72      1
#> [283,]      39     64
#> [284,]      72     54
#> [285,]      65     57
#> [286,]     347     31
#> [287,]     312     50
#> [288,]      93     17
#> [289,]      61     25
#> [290,]      65     28
#> [291,]      54      1
#> [292,]      51      5
#> [293,]     252      5
#> [294,]     255     10
#> [295,]      86     20
#> [296,]      81     22
#> [297,]      64     13
#> [298,]     245     11
#> [299,]     244      8
#> [300,]      65      6
#> [301,]      63      3
#> [302,]      67      1
#> [303,]      31     20
#> [304,]      73     56
#> [305,]      49     40
#> [306,]     111     25
#> [307,]      66     22
#> [308,]     293     10
#> [309,]      71     20
#> [310,]      51     29
#> [311,]      36     45
#> [312,]      70     25
#> [313,]      67      7
#> [314,]      60     20
#> [315,]      53     16
#> [316,]      71     12
#> [317,]      74     21
#> [318,]     261     17
#> [319,]      49     30
#> [320,]      61     13
#> [321,]      57     39
#> [322,]     218      2
#> [323,]      54     15
#> [324,]     111     53
#> [325,]      76     14
#> [326,]      78      0
#> [327,]      77     48
#> [328,]      83      8
#> [329,]     253      7
#> [330,]     254     12
#> [331,]      72      3
#> [332,]      71      4
#> [333,]     139      2
#> [334,]      61     18
#> [335,]      76      4
#> [336,]     125     22
#> [337,]     124     26
#> [338,]     194     53
#> [339,]     225     13
#> [340,]      62      8
#> [341,]     230      2
#> [342,]      78      9
#> [343,]     237     12
#> [344,]      44      1
#> [345,]     248     15
#> [346,]     258      4
#> [347,]     232     14
#> 
```
