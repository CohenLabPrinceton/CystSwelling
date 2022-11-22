time_x = 0;
time_y = 25;
time_font = 40;
field_x = time_x;
field_y = time_y + time_font + time_font/2.5;
field_font = 40;

//setForegroundColor(3, 3, 3); // if black
setForegroundColor(251, 251, 251); // if white
run("Time Stamper", "starting=0 interval=10 x="+toString(time_x)+" y="+toString(time_y)+" font="+time_font+" decimal=0 anti-aliased or=min"); 
// xy 50 50 font 150 for 2048x2048 image
run("Scale Bar...", "width=50 height=8 font=30 color=Black background=None location=[Lower Right] bold overlay");
 	
for (i = 1; i <= 6; i++) {
	setSlice(i);
	setFont("SansSerif", field_font, " antialiased");
	setColor("#FFFFFF"); // for black: 000000, white: FFFFFF
	drawString("Field OFF", field_x, field_y); // x, y of bottom left corner
	//run("Draw", "slice");
}
// font 100 xy 100 200 for  2048x2048 image

for (i = 7; i <= 49; i++) {
	setSlice(i);
	setFont("SansSerif", field_font, " antialiased");
	setColor("#FFFFFF");
	drawString("Field -->", field_x, field_y); // x, y should be same as x from time stamper; y is y of timestamper + font size
	//run("Draw", "slice");
}

