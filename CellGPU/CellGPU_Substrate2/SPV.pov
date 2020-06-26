#include "colors.inc"
#include "finish.inc"
#include "shapes.inc"
#include "textures.inc"

global_settings {
   adc_bailout 0.0039216
   ambient_light rgb <1,1,1>
   assumed_gamma 1.5
   noise_generator 2
}

#macro Colloid2(Pos)
sphere{Pos, 0.5 texture{pigment{ color rgb <0.2,0.2,0.8> transmit 0.25} finish {ambient 0.1 reflection 0.0 phong 0.25}}}
#end //-- end of macro Colloid2(Pos)

#macro Voro(Pos1, Pos2, Pos3, Pos4, Pos5, Pos6, Color1, Color2, Color3)

triangle{ <Pos1, Pos2, 0>, <Pos3, Pos4, 0>, <Pos5, Pos6, 0> texture{pigment{color rgb <Color1,Color2,Color3>} finish {ambient 0.1 reflection 0.0 phong 0.1}}}
cylinder{<Pos1, Pos2, -0.1>, <Pos5,Pos6,0.1>, 0.01 texture{pigment{ color rgb <0,0,0> transmit 0.0} finish {ambient 0.1 reflection 0.0 phong 0.1}}}

#end

#macro VoroC(Pos1, Pos2, R, Color1, Color2, Color3)

//cylinder{<Pos1, Pos2, -0.1>, <Pos1,Pos2, -0.000001>, R texture{pigment{ color rgb <Color1, Color2, Color3> transmit 0.1} finish {ambient 0.1 reflection 0.0 phong 0.1}}}
cylinder{<Pos1, Pos2, -0.1>, <Pos1,Pos2, 0.1>, 0.1 texture{pigment{ color rgb <1,1,1> transmit 0.0} finish {ambient 0.1 reflection 0.0 phong 0.1}}}

#end

#macro VoroL(Pos1, Pos2, Pos3, Pos4, Pos5, Pos6, Color1, Color2, Color3, R)

intersection
{
triangle{ <Pos1, Pos2, 0>, <Pos3, Pos4, 0>, <Pos5, Pos6, 0> texture{pigment{color rgb <Color1,Color2,Color3> } finish {ambient 0.1 reflection 0.0 phong 0.1}}}
//prism{-0.5, 0.5, 4 <Pos1, Pos2>, <Pos3, Pos4>, <Pos5, Pos6>, <Pos1, Pos2> texture{pigment{color rgb <Color1,Color2,Color3> } finish {ambient 0.1 reflection 0.0 phong 0.1}} rotate <-90,0,0> }
cylinder{ <Pos3, Pos4, 0.1>, <Pos3,Pos4,-0.1>, R texture{pigment{ color rgb <Color1,Color2,Color3>} finish {ambient 0.1 reflection 0.0 phong 0.1}}}
}

#end

#macro VoroD(Pos1, Pos2, Pos3, Pos4, Pos5, Pos6, Color1, Color2, Color3, R)

difference
{
cylinder{ <Pos3, Pos4, 0.1>, <Pos3,Pos4,-0.1>, R texture{pigment{ color rgb <Color1,Color2,Color3>} finish {ambient 0.1 reflection 0.0 phong 0.1}}}
prism{-0.5, 0.5, 4 <Pos1, Pos2>, <Pos3, Pos4>, <Pos5, Pos6>, <Pos1, Pos2> texture{pigment{color rgb <Color1,Color2,Color3> } finish {ambient 0.1 reflection 0.0 phong 0.1}} rotate <-90,0,0> }
//triangle{ <Pos1, Pos2, 0>, <Pos3, Pos4, 0>, <Pos5, Pos6, 0> texture{pigment{color rgb <Color1,Color2,Color3> } finish {ambient 0.1 reflection 0.0 phong 0.1}}}
}

#end

#macro VoroA(Pos1, Pos2, R, Color1, Color2, Color3)

cylinder{<Pos1, Pos2, -0.1>, <Pos1,Pos2, 0.1>, 0.1 texture{pigment{ color rgb <1,1,1> transmit 0.0} finish {ambient 0.1 reflection 0.0 phong 0.1}}}
cylinder{<Pos1, Pos2, -0.1>, <Pos1,Pos2, -0.000001>, R texture{pigment{ color rgb <Color1, Color2, Color3> transmit 0.1} finish {ambient 0.1 reflection 0.0 phong 0.1}}}

#end


//Licht
light_source { <-200,150,-200>  color Grey shadowless }
light_source { <-200,150,200>  color Grey shadowless }
light_source { <200,150,-200>  color Grey shadowless }
light_source { <200,150,200>  color Grey shadowless }

# declare campos = <5,5, 12>;

# declare camlook = <5,5, 0>;

background {
  //color rgb <0., 0., 0. >
  color rgb <1,1,1 >
}

/*
light_source {
  <128,-128, -128> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <128, 128, 128> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <0, 128, 128> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <128, 0, 128> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <128, 128, 0> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <0, 0, 128> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <0, 128, 0> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <128, 0, 0> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <0, 0, 0> 
  color rgb <1,1.0,1>
  spotlight
  radius 100
  falloff 40
  tightness 30
  point_at <64, 64, 64>
}
light_source {
  <0, -550, 900> 
  color rgb <1,1.0,1>
  spotlight
  radius 50
  falloff 40
  tightness 30
  point_at <0.0, 100.0, 250.0>
}
light_source {
  <0, -550, 2000> 
  color rgb <1,1.0,1>
  spotlight
  radius 50
  falloff 40
  tightness 30
  point_at <0.0, 00.0, 550.0>
}
light_source {
  <0, -600, 200> 
  color rgb <1,1.0,1>
  spotlight
  radius 50
  falloff 40
  tightness 30
  point_at <0.0, 0.0, 200.0>
}

light_source {
  <100, 600, -100> 
  color rgb <0.8,1.0,0.8>
  spotlight
  radius 50
  falloff 40
  tightness 30
  point_at <0.0, 300.0, -0.0>
}

light_source {
  <0, 600, 200>
  color rgb  <0.8,1.0,0.8>
  spotlight
  radius 50
  falloff 40
  tightness 10
  point_at <0.0, 300, 0.0>
  shadowless
}
light_source {
  <0, 600, 0>
  color rgb  <0.8,1.0,0.8>
  spotlight
  radius 50
  falloff 40
  tightness 10
  point_at <0.0, 0,-50.0>
  shadowless
}
light_source {
  <0, 600, 0>
  color rgb <0.8,1.0,0.8>
  spotlight
  radius 50
  falloff 40
  tightness 10
  point_at <0.0, 0, -100.0>
  shadowless
}

light_source {
  <0, 10, 400>
  color rgb <0.5,0.6,0.5>
  spotlight
  radius 50
  falloff 40
  tightness 10
  point_at <0.0, 200, 0.0>
  shadowless
}

light_source {
  <100, -100, 153>
  color rgb <0.5,0.6,0.5>
  spotlight
  radius 1000
  falloff 40
  tightness 10
  point_at <0.0, 0, 0.0>
  shadowless
}
*/

camera {
location  campos
look_at camlook
up    <0,1,0>
right  <1,0,0>
}

/*
#switch(frame_number)
#range(1, 9)
#include concat("snap00", str(frame_number, -1, 0), ".pov")
#break
#range(10, 99)
#include concat("snap00", str(frame_number, -2, 0), ".pov")
#break
#range(100, 999)
#include concat("snap00", str(frame_number, -3, 0), ".pov")
#break
#range(1000, 9999)
#include concat("snap00", str(frame_number, -4, 0), ".pov")
#break
#end
*/

#include "snap001.pov"

