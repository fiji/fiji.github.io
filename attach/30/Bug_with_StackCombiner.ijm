newImage("MyRamp", "8-bit Ramp", 600, 50, 1);
run("Rotate 90 Degrees Right");

newImage("MyImage", "8-bit Ramp", 600, 200, 1);
run("Flip Horizontally");
run("Rotate 90 Degrees Right");

run("Stack Combiner", "stack=MyImage stack=MyRamp");