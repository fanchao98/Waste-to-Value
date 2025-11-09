# Waste-to-Value: Reutilized Material Maximization for Additive and Subtractive Hybrid Remanufacturing [[Project page]](https://fanchao98.github.io/Refab%20page/refab.html)
![.](/teaser.png)
Turning damaged or discarded parts into high-value components unlocks hidden potential and gives waste a second life with purpose. Starting with the initial worn Bearing model (green) and the target topologically optimized Bearing model (blue) (a), our computational framework generates process planning solution aimed at maximizing material savings by reutilizing as much volume as possible while ensuring manufacturability. Excess material is removed via subtractive manufacturing (b), and additional material is subsequently added through additive manufacturing (c) to obtain the target model. Compared to fabricating the target model entirely using additive manufacturing, our remanufacturing process planning solution achieve material savings of up to 90%.
# Dependency
You should install **CGAL(5.5.1)**, **opencv**, **libigl**, **boost(1_80_0)**, and include glpk-5.0\src, before compiling this project. I recommend using Visual Studio 2022.
# Start
You can simply place the input initial model (_ori.off) and target model (_target.off) in the "models" folder, and run.

# License
All rights about the program are reserved by the authors of this project. The programs can only be used for research purpose. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program.
