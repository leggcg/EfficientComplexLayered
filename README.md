# An Efficient Transport Estimator for Complex Layered Materials

[Luis E. Gamboa](https://lc.fie.umich.mx/~legg/), [Adrien Gruson](https://beltegeuse.github.io/research/), [Derek Nowrouzezahrai](http://www.cim.mcgill.ca/~derek/).
In Computer Graphics Forum (Eurographics 2020).
[[Project page]](https://lc.fie.umich.mx/~legg/complexlayered.php)

## Installation
This is a branch of mitsuba 0.6.0 renderer.
In this release, we provide the following plugins inside `src/devplugins`:
- `fast_accurate`: our method
- `multilayered_guo2018`: Guo2018 method ([project page](https://shuangz.com/projects/layered-sa18/)) 
- and all the required supporting plugins to render with both methods.


## Scenes
- <a href="https://lc.fie.umich.mx/~legg/papers/complexlayered/scenes/coffeetable.zip">Coffee table, including the red mug.</a>

## Rendering
To render a scene file using our material set and 512 samples per pixel use:
`mitsuba -D materials=ours -D spp=512 scene.xml`

Set `materials=guo` to render using Guo2018.

To test different integrators, use `-D integrator=` to:
- `pathwsampler`: path tracing w/MIS.
- `directwsampler`: direct illumination w/MIS.
- `directlightsampling`: direct illumination using only Light Sampling.
