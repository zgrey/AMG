# Python Program for quick gmsh grid gen

import string,sys

def convert_pts_togmsh(filename,outputname,startpoint):
	
	# Read in data using this bit
	fin=open(filename, 'r')
	i=0
	x = []; y = []; z = [];
	
	lines = fin.readlines()

	for line in lines:
		# comments indicated by #
		if line[0] != '#' or line[0] != '':
			i=i+1
			data = string.split(line)
			x.append(float(data[0]))
			y.append(float(data[1]))
			z.append(float(data[2]))

			n_lines = int(i)

	# Write data out with this;

	fout = open(outputname,'w')

	lc_name = "%s_lc" % filename[0:3]
	# Format
	# Point(1) = {0, 0, 0, lc};
	fout.write("%s = 0.0075;\n" % lc_name)
	j = startpoint
	for i in range(n_lines):
		outputline  = "Point(%i) = { %8.8f, %8.8f, %8.8f, %s};\n " \
					      % (j,x[i],y[i],z[i],lc_name )
		j = j + 1
		fout.write(outputline)

	# gmsh bspline format	
	# Write out splinefit line
	# Suction Side
	fout.write("BSpline(%i) = {%i:%i};\n" \
		   % (startpoint,startpoint,startpoint+(n_lines-1)*0.5))
		   
	fout.write("BSpline(%i) = {%i:%i,%i};\n" \
		   % (startpoint+1,startpoint+(n_lines-1)*0.5,startpoint+n_lines-2,startpoint))
	
	# Write out Boundaries
	fout.write("radius = 4;\n")
	# used for first set of meshes
	# fout.write("bc_lc = 0.25;\n")
	fout.write("bc_lc = 0.1;\n")
	fout.write("Point(1) = {0.5, 0, 0, bc_lc};\n")
	fout.write("Point(2) = {radius+0.5, 0, 0, bc_lc};\n")
	fout.write("Point(3) = {-radius+0.5, 0, 0, bc_lc};\n")
	#fout.write("Point(4) = {0.5, -radius, 0, bc_lc};\n")
	#fout.write("Point(5) = {0.5, radius, 0, bc_lc};\n")
	
	#fout.write("Circle(1) = {5, 1, 2};\n") 
	#fout.write("Circle(2) = {2, 1, 4};\n")
	#fout.write("Circle(3) = {4, 1, 3};\n")
	#fout.write("Circle(4) = {3, 1, 5};\n")
	
	fout.write("Circle(1) = {2, 1, 3};\n") 
	fout.write("Circle(2) = {3, 1, 2};\n")
	
	#fout.write("Line Loop(5) = {1, 2, 3, 4};\n")
	fout.write("Line Loop(5) = {1, 2};\n")
	fout.write("Line Loop(6) = {%i,%i};\n" \
		   % (startpoint,startpoint+1))
		   
	# Construct Surface
	fout.write("Plane Surface(%i) = {%i,%i};\n" \
		   % (startpoint+1,5,6))
	#fout.write("Physical Line('FARFIELD') = {1, 2, 3, 4};\n")
	fout.write("Physical Line('FARFIELD') = {1, 2};\n")
	fout.write("Physical Line('airfoil') = {%i,%i};\n" \
		   % (startpoint,startpoint+1))
	fout.write("Physical Surface('fluid') = {%i};\n" \
		   % (startpoint + 1))
	
	# Close files
	fout.close
	fin.close
	
def main():
  inputfile = sys.argv[1]
  convert_pts_togmsh(inputfile,inputfile+".geo",1000)

if __name__ == "__main__":
    main()
