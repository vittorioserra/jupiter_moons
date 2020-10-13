#include "colors.inc"
#include "textures.inc"



camera
{
        location <-3000000.0, 500000.0, -1000000.0>
        look_at <0.0, 0.0,  0.0>
        right x*image_width/image_height
}
#declare sky_scale = 0.002;

//SUN unique light source
light_source
{
   	 <-778510000, 0 , 0>
   	 color White
   	 looks_like { sphere{<0,0,0> 1392700}}
}

//ALL dimensions of the planets are reported in km

// Jupiter
sphere
{
        <0,0,0> 142984
        texture { pigment{ image_map { png "jupiter_texture.png" map_type 1 } } }
        rotate<0,-90,0>
}

//Moons of Jupiter

//Io
sphere
{
        <421700,0,0> 3643.2
        texture { pigment{ image_map { png "io_surface_texture.png" map_type 1 } } }
        rotate<0,-90,0>
}

//Europa
sphere
{
        <0,0,671034> 3121.6
        texture { pigment{ image_map { png "europa_surface_texture.png" map_type 1 } } }
        rotate<0,-90,0>
}

//Ganymede
sphere
{
        <1070412,0,0> 5262.4
        texture { pigment{ image_map { png "ganymede_surface_texture.png" map_type 1 } } }
        rotate<0,-90,0>
}

//Callisto
sphere
{
        <0,0,1882709> 4820.6
        texture { 
		
		pigment{
			//rgb<1,0,0>
			image_map { png "callisto_surface_texture.png" map_type 1 } 
			} 
	}
        rotate<0,-90,0>
}


// Sky
sky_sphere
{
        pigment
        {
                crackle form <1,1,0> color_map { [.3 rgb 1] [.4 rgb 0] }
                scale sky_scale
        }
}







