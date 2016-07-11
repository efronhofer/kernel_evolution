program kernelEvolution;

// Copyright (C) 2015  Emanuel A. Fronhofer
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// additional libraries
uses Math;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// CONSTANTS
const
	X_MAX = 100;														// world dimensions
	Y_MAX = 100;
	CELL_MAX = 100;														// max no. of trees per cell
	KERNEL_MAX = 21;													// maximum kernel distance

	PI = 3.141593;														// pi
	
	DELTA_MUTATION = 2;													// mutation width
	SEG_SD = 0.1;														// segregation standard deviation for segregation kernel
	COMPETITION_THRESHOLD = 0.01;										// threshold for calculation of competition
	INIT_TCELL = 1;														// initial no of trees per cell
	
	//RS = 2;															// random seed
	
	TIME_MAX = 10000;													// simulation time
	RUN_MAX = 25;														// no of repeats
	
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// TYPE DECLARATION
type
	//record for a single tree
	TTree = record
		age : word;
		kernel : array[1..KERNEL_MAX] of single;						// phenotypic kernel
		location : array[1..2] of single;
		mort : real;
	end;
	
	//record for a cell
	TCell = record
		tree : array[1..CELL_MAX] of TTree;								// trees in one cell
		no_trees : word;												// present no of trees
	end;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// GLOBAL VARIABLES
var
	cell : array[1..X_MAX, 1..Y_MAX] of TCell;							// all trees in the world
	time : longint;														// counter for main time loop
	
	fertility : real;													// fertility per individual and time-step
	mort_0 : real;														// baseline mortality (density independent)
	hsc_age : real;														// half saturation constant of age related competition
	mort_disp0 : real;													// per step dispersal mortality
	mutation_rate : real;												// mutation rate for kernel
	comp_kernel_const : real;											// competition kernel constant
	sigma : real;														// standard deviation of competition kernel
	kurtosis : real;													// kurtosis of competition kernel
	
	population_size : longint;											// population size (for per generation output)
	population_size_afterComp : longint;								// population size (for per generation output) after competition
	mean_dispDist : real;												// mean dispersal distance (for per generation output)
	mean_dispDist_cnt : longint;										// counter for no indiv used to calc mean disp dist
	compDist : integer;													// maximal distance at which competition is calculated (cells)
	run_no : integer;													// actual repeat
	
	dispCostCheckBox, dispControlCheckBox: string;						// determines whether dispersal costs are kernel (maternal, i.e. trade-off) or distance dependent (offspring)
																		// determines whether dispersal is determiend by mother's or offspring kernel
	
	directions : array[1..2,1..4] of integer;							// help variable for mating with nearest neighbour; search procedure
	
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// FUNCTIONS

// ---------------------------------------------------------------------
// log Gauss distribution ----------------------------------------------
function loggauss(const mw, sd: real): real;
// calculates the mean random fertility
// of a single female individual from a
// Gaussian probability distribution,
// using the variable "Fertility" as
// mean value and "sigma" (environmental
// fluctuations) as standard deviation. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
	end;

//----------------------------------------------------------------------
// gamma function ------------------------------------------------------
function gamma(x: real): real;
	// gamma function
	// From Borland Pascal Programs for Scientists and Engineers
	// by Alan R. Miller, Copyright C 1993, SYBEX Inc
	end;

// --------------------------------------------------------------------------
// poisson distribution -----------------------------------------------------
function Poisson(lamb:double):integer;
// this function creates Poisson distributed random numbers
// with mean lamb. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
	end;

//--------------------------------------------------------------------------
// gaussian normal distribution ---------------------------------------------
function gauss : real;
// this function creates Gaussian distributed random numbers with mean
// 0 and standard deviation 1. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
	end;

// --------------------------------------------------------------------------
// mutation -----------------------------------------------------------------
function inheritance(const allele : real) : real;
	// phenotypic model; i.e. mid-parent value and segregation kernel (standard deviation of gaussian dist. SEG_SD)
	// introduces additionally a mutation with a gaussian distribution around the mean "allele" with
	// standard deviation DELTA_MUTATION (see constants)!
	// called by procedure "dispersal"
	var
		mut_sd : real;			// delta mutation as a function of time
	
	begin
		// calculate additional mutation sd; decreases over time; stabilizing selection
		mut_sd := 1;
		// mutation or not 
		if random < mutation_rate then
		begin
			mut_sd := DELTA_MUTATION * exp(-5*time/TIME_MAX);
			inheritance := allele * loggauss(1,SEG_SD+mut_sd);
		end
		else
		begin
			inheritance := allele * loggauss(1,SEG_SD);
		end;
	end;

//----------------------------------------------------------------------
// checks coordinates and makes a torus --------------------------------
function maketorus(const x: integer):integer;
begin
	maketorus := x;
	
	if x > X_MAX then
	begin
		maketorus := x mod X_MAX;
		if x = 0 then maketorus := X_MAX;
	end;

	if x < 1 then
	begin
		maketorus := X_MAX - (abs(x) mod X_MAX);
	end;
end;

// ---------------------------------------------------------------------
// the competition kernel ----------------------------------------------
function comp_kernel (const d : real): real;
	// as in Roughgarden 1974 
	// comp_kernel_const is calculated once at simulation start in the initialization procedure!
	begin
		comp_kernel := exp(-power((abs(d)/comp_kernel_const),kurtosis));
	end;

//----------------------------------------------------------------------
// calculate mortality by competition ----------------------------------
function mort_comp(const x_act,y_act, t_act : word): real;
var
	x,y : integer;														// counter for x and y dimensions
	xl, yl : integer;
	t : word;															// counter for trees in cell
	xcmin, xcmax, ycmin, ycmax : integer;								// area of competition
	prod_surv : real;													// product of 1-comp i.e. survival probabilities
	indiv_comp : real;													// competition between actual and focal tree
	comp_dist : real;													// distance between actual and focal tree
	
begin
	// initialize sum_comp
	prod_surv := 1;
	// calculate area of competition (cells)
	xcmin := x_act - (compDist-1);
	xcmax := x_act + (compDist-1);
	ycmin := y_act - (compDist-1);
	ycmax := y_act + (compDist-1);
	// sum up all competitive interactions with tree t in cell x,y; only if the focal tree is competitive, i.e. was alive at the beginning
	// of this time step; trees are competitive as long as the actual time step is running, even if they die
	// this avoids artefacts of tree order!
	for x := xcmin to xcmax do
	for y := ycmin to ycmax do
	begin
		xl := maketorus(x);
		yl := maketorus(y);
		for t := 1 to cell[xl,yl].no_trees do
		begin
			if not((xl = x_act)and(yl = y_act)and(t = t_act)) then
			begin
				//distance between focal and actual tree 
				comp_dist := sqrt((cell[xl,yl].tree[t].location[1]-cell[x_act,y_act].tree[t_act].location[1])*(cell[xl,yl].tree[t].location[1]-cell[x_act,y_act].tree[t_act].location[1]) + (cell[xl,yl].tree[t].location[2]-cell[x_act,y_act].tree[t_act].location[2])*(cell[xl,yl].tree[t].location[2]-cell[x_act,y_act].tree[t_act].location[2]));
				// calculate individual part of competition (age and distance)
				indiv_comp := (cell[xl,yl].tree[t].age / (cell[xl,yl].tree[t].age + hsc_age) ) * comp_kernel(comp_dist);
				// add to sum
				prod_surv := prod_surv * (1-indiv_comp);
			end;
		end;
	end;
	
	mort_comp :=  1-prod_surv;
end;

//----------------------------------------------------------------------
// dispersal distance --------------------------------------------------
function distance (const kernel: array of single) : real;
var
	k : integer;														// counter for kernel distance classes
	threshold : array[0..KERNEL_MAX] of single;							// calc. thresholds for distance
	norm_fact : single;													// normalization factor for kernel
	limit : single;														// random no. for dispersal distance
	int_dist : integer;													// integer dispersal distance taken from kernel
	mean_dispDist_act : real;											// mean dispersal distance of the focal animal
	sum_kernelProbs : real;												// sum of ketnel probabilities
	helpdist : real;													// help dist var for area to area disp
	
begin
	// initialize focal mean dispersal distance and sum
	mean_dispDist_act := 0;
	sum_kernelProbs := 0;
	// initialize threshold array
	threshold[0] := 0;
	// note that the counter has to start with 0 although the distance class is one!
	for k:=0 to (KERNEL_MAX-1) do threshold[k+1] := threshold[k]+kernel[k];
	// now normalize threshold between 0 and 1, i.e. set last value to 1 and recalc the others
	norm_fact := threshold[KERNEL_MAX];
	if norm_fact = 0 then
	begin
		writeln('ERROR: no kernel!');
	end
	else
	begin
		
		for k:=1 to (KERNEL_MAX) do 
		begin
			threshold[k] := threshold[k]/norm_fact;
			// calculate mean dispersal distance for output
			mean_dispDist_act := mean_dispDist_act + (k-1)*kernel[k-1];
			// sum up kernel probabilities
			sum_kernelProbs := sum_kernelProbs + kernel[k-1];
		end;
		
		// calculate mean dispersal distance for output
		mean_dispDist_act := mean_dispDist_act/sum_kernelProbs;
		mean_dispDist := mean_dispDist + mean_dispDist_act;
		inc(mean_dispDist_cnt);
		
		// draw a random number
		limit := random;
		// determine dispersal distance
		for k:=1 to KERNEL_MAX do
		begin
			if limit <= threshold[k] then
			begin
				int_dist := k;
				break;
			end;
		end;
		
		//include cat zero
		int_dist := int_dist -1;
		
		// this ony returns integer numbers; thus I include area-to-area dispersal here
		if int_dist = 0 then helpdist := 0
		else helpdist := int_dist-random;
		//helpdist := int_dist+random;
		
		if helpdist < 0 then writeln('ERROR: dispersal distance under zero!');
		
		distance := helpdist;
	end;
end;

//----------------------------------------------------------------------
// distance dependent dispersal mortality ------------------------------
// mortality probability after x steps follows a geometric distribution
// in addition we assume a random walk (kernel), i.e. the really covered 
// distance is proportional to the net distance squared
function mort_disp (const d:real): real;
begin
	mort_disp := 1 - exp(-mort_disp0*d);
end;

//----------------------------------------------------------------------
// kernel dependent dispersal mortality --------------------------------
function mort_kernel (const kernel: array of single):real;
var
	k : integer;														// counter for kernel
	zaehler, nenner : real;												// hilfsvariablen für zähler und nenner
begin
	zaehler := 0;
	nenner := 0;
	for k := 0 to (KERNEL_MAX-1) do 
	begin
		if k = 0 then zaehler := zaehler + kernel[k] * mort_disp(k)
		else zaehler := zaehler + kernel[k] * mort_disp(k-0.5);
		nenner := nenner + kernel[k];
	end;
	mort_kernel := zaehler/nenner;
end;



//----------------------------------------------------------------------
//----------------------------------------------------------------------
// PROCEDURES

//----------------------------------------------------------------------
// initialization ------------------------------------------------------
procedure initialize;
var
	x, y: word;															// counter for x and y dimensions
	k : word;															// counter for kernel initialization
	t : word;															// counter for trees in a cell
	
	input : text; 														// data input file

begin
	// initialize random number generator
	//randseed := RS;
	randomize;
	
	// parameter einlesen
	assign(input, 'input_kernelEvolution.in');
	reset(input);
		
	readln(input); readln(input); readln(input);						// skip header
	readln(input, fertility); readln(input);							// read fertility
	readln(input, mort_0); readln(input);								// read morality
	readln(input, hsc_age);readln(input);								// half saturation constant for age related competition
	readln(input, sigma); readln(input);								// standard deviation of distance related competition kernel
	readln(input, kurtosis); readln(input);								// kurtosis of distance related competition kernel
	readln(input, mort_disp0); readln(input);							// per step dispersal mortality
	readln(input, mutation_rate); readln(input);						// mutation rate for kernel
	readln(input, dispCostCheckBox); readln(input);						// determines whether dispersal costs are kernel or distance dependent
	readln(input, dispControlCheckBox); readln(input);					// determines whether dispersal is controlled by mother or offspring
		
	close(input);
	
	// initialize trees in world
	for x:= 1 to X_MAX do
	for y:= 1 to Y_MAX do
	begin // world loop
		cell[x,y].no_trees := 0;
		
		for t:=1 to INIT_TCELL do
		begin
			// all trees are initialized "alive"
			cell[x,y].tree[t].mort := 0;
		
			// initialize age
			cell[x,y].tree[t].age := 1;
			
			//initialize location
			cell[x,y].tree[t].location[1] := x-1 + random;
			cell[x,y].tree[t].location[2] := y-1 + random;
		
			// initialize kernel
			for k := 1 to KERNEL_MAX do 
			begin
				cell[x,y].tree[t].kernel[k] := 0.5 + random*0.2-0.1;
			end;
			
			// increase tree counter
			inc(cell[x,y].no_trees);
		end; // cell loop
	end; // world loop
	
	// calculate constant for competition kernel
	comp_kernel_const := (sigma*sqrt(gamma(1/kurtosis))) / sqrt(gamma(3/kurtosis));
	// calculate no of cells used for competition
	compDist:= ceil( power((-ln(COMPETITION_THRESHOLD)),(1/kurtosis))*(sigma*sqrt(gamma(1/kurtosis))) / sqrt(gamma(3/kurtosis)) );
	
	// init directions
	directions[1,1] := 0;
	directions[1,2] := -1;
	directions[1,3] := 0;
	directions[1,4] := 1;
	directions[2,1] := -1;
	directions[2,2] := 0;
	directions[2,3] := 1;
	directions[2,4] := 0;

end;

//----------------------------------------------------------------------
// analysis and data output --------------------------------------------
procedure analysis;
var
	x,y : word;															// counter for x and y dimensions
	t : word;															// counter for trees
	k : word;															// counter for kernel output
	
	output : text;														// data output file
	outputfilename, runstring : string;									// helpstings for output file path

begin
	// create output file
	str(run_no, runstring);
	outputfilename := './results/output_kernelEvolution_run' + runstring + '.out'; 
	assign(output, outputfilename);
	rewrite(output);
	
	// write output header
	write(output, '     x		locx		   y		locy		age');
	for k := 1 to KERNEL_MAX do write(output,'            k',k);
	writeln(output);

	// loop over all trees
	for x := 1 to X_MAX do
	for y := 1 to Y_MAX do
	begin //1
		for t := 1 to cell[x,y].no_trees do
		begin //2
			
				write(output, x:6, '    ',cell[x,y].tree[t].location[1]:10:2,'     ', y:6, '    ',cell[x,y].tree[t].location[2]:10:2,'     ', cell[x,y].tree[t].age:6);
				for k := 1 to KERNEL_MAX do write(output,'     ', cell[x,y].tree[t].kernel[k]:10:5);
				writeln(output);
	
		end; //2
	end; //1
	
	// close output file
	close(output);
end;

//----------------------------------------------------------------------
// death procedure -----------------------------------------------------
// effects of age and competition
procedure death;
var
	x, y : word;														// counter for x and y dimensions
	t : word;															// counter for trees in a cell
	
begin
	// reset counter for population size
	population_size := 0;
	
	for x := 1 to X_MAX do
	for y := 1 to Y_MAX do
	begin //1
		// density regulation		
		for t := 1 to cell[x,y].no_trees do
		begin //2
			// increase counter for population size
			inc(population_size);
			cell[x,y].tree[t].mort := 1 - (1 - mort_0) * (1 - mort_comp(x,y,t) );
		end; //2
	end; // 1
	
	// now really kill a dead individuals
	population_size_afterComp :=0;
	for x:=1 to X_MAX do
	for y:=1 to Y_MAX do
	begin
		t := 0;
		while t < cell[x,y].no_trees do
		begin //2
			inc(t);
			// apply mortality and increase age
			if random < cell[x,y].tree[t].mort then
			begin
				cell[x,y].tree[t] := cell[x,y].tree[cell[x,y].no_trees];
				dec(cell[x,y].no_trees);
				dec(t);
			end
			else 
			begin
				inc(cell[x,y].tree[t].age);
				inc(population_size_afterComp);	
			end; //3
			
		end; //2
	end;
	
end;

//----------------------------------------------------------------------
// dispersal procedure called by reproduction procedure ----------------
procedure dispersal(const x,y,t,xf,yf,tf : integer);
var
	disp_dist : real;													// dispersal distance of focal offspring
	phi : real;															// angle for dispersal
	locendx, locendy : real;											// end location (plus random for grid)
	xend, yend : integer;												// end coordinates adapted to grid
	diffx, diffy : real;												// difference between cell and location
	disp_costs : real;													// dispersal costs
	k : word;															// counter for loop over kernel
	helpkernel : array[1..KERNEL_MAX] of single;						// help variable for offspring kernel
	
begin
	// let's start with inheritance (note: phenotypic model; i.e. mid-parent value and segregation kernel):
	for k:=1 to KERNEL_MAX do
	begin
		helpkernel[k] := inheritance((cell[x,y].tree[t].kernel[k]+cell[xf,yf].tree[tf].kernel[k])/2);
	end;

	// determine new location for offspring from dispersal kernel of the mother tree
	// first determine distance
	if (dispControlCheckBox = 'maternal') then disp_dist := distance(cell[x,y].tree[t].kernel);
	if (dispControlCheckBox = 'offspring') then disp_dist := distance(helpkernel);
	if (not(dispControlCheckBox = 'maternal') and not(dispControlCheckBox = 'offspring')) then writeln('ERROR: dispControlCheckBox must be one of mother or offspring');
	
	// now determine angle; randomly drawn for 2pi
	phi := 2 * PI * random;
	// location end
	locendx := cell[x,y].tree[t].location[1] + disp_dist * cos(phi);    
	locendy := cell[x,y].tree[t].location[2] + disp_dist * sin(phi);
	// target coordinates
	xend := ceil( locendx );    
	yend := ceil( locendy );
	diffx := locendx-xend;
	diffy := locendy-yend;
	
	// make torus
	xend := maketorus(xend);
	yend := maketorus(yend);
	locendx := xend  + diffx;
	locendy := yend  + diffy;

	// check whether tree is empty
	if cell[xend,yend].no_trees <= CELL_MAX then
	begin //1
		// take into account dispersal costs; these may either be distance or kernel dependent
		if dispCostCheckBox = 'offspring' then disp_costs := mort_disp(disp_dist)
		else
		begin
			if dispCostCheckBox = 'maternal' then 
			begin
				disp_costs := 0;
				//if (dispControlCheckBox = 'mother') then disp_costs := mort_kernel(cell[x,y].tree[t].kernel);
				//if (dispControlCheckBox = 'offspring') then disp_costs := mort_kernel(helpkernel);
			end
			else writeln('ERROR: dispCostCheckBox must be one of KERNEL or DISTANCE!');
		end;
		
		// check whether the individual survives the dispersal event
		if random >= disp_costs then
		begin //2
			// inc no of trees in target cell
			inc(cell[xend,yend].no_trees);
			
			// here include inheritance etc.
			cell[xend,yend].tree[cell[xend,yend].no_trees].age := 1;

			cell[xend,yend].tree[cell[xend,yend].no_trees].kernel := helpkernel;
			
			cell[xend,yend].tree[cell[xend,yend].no_trees].location[1] := locendx;
			cell[xend,yend].tree[cell[xend,yend].no_trees].location[2]:= locendy;
		end; //2
	end //1
	else writeln('ERROR: cell is overfull!');
			
end;

//----------------------------------------------------------------------
// reproduction procedure ----------------------------------------------
// includes reproduction and dispersal
procedure reproduction;
var
	x,y : word;															// counter for x and y dimensions
	t : word;															// counter for tree t in cell
	xf,yf,tf : word;													// address of the "male"
	o : word;															// counter for no. of offspring
	no_offspring : word;												// no. of offspring of the focal tree
	stop : boolean;
	potmates : array[1..3,1..1000] of integer;							// help vect for animals
	potmates_dist : array[1..1000] of real;								// help vect for animals distance
	potmates_cnt,i,j,l,m, pm, nearestneighbour, ring : word;			// counter for potential mates
	mindist : real;														// help var to save minimum distance for mate choice
	xr,yr,dir,hxr,hyr : integer; 										//for ring search
	tradeOff_costs: real;												// trade off

begin
	for x:= 1 to X_MAX do
	for y:= 1 to Y_MAX do
	begin
		for t := 1 to cell[x,y].no_trees do
		begin //2
			// reproduction only for adults!
			if cell[x,y].tree[t].age > 1 then
			begin //3
				// mate locally: check the nn8 etc rings and then choose randomly
				stop := false;
				ring := 0;
				potmates_cnt := 0;
				// first check whether there is another tree in the present cell
				if cell[x,y].no_trees > 1 then
				begin //4
					// do not mate with new indivs from this time step
					for j:= 1 to cell[x,y].no_trees do
					begin //5
						if not(j=t) then
						begin //6
							if cell[x,y].tree[j].age > 1 then
							begin
							inc(potmates_cnt);
							potmates[1,potmates_cnt] := x;
							potmates[2,potmates_cnt] := y;
							potmates[3,potmates_cnt] := j;
							potmates_dist[potmates_cnt] := (cell[x,y].tree[t].location[1] - cell[x,y].tree[j].location[1])*(cell[x,y].tree[t].location[1] - cell[x,y].tree[j].location[1]) + (cell[x,y].tree[t].location[2] - cell[x,y].tree[j].location[2])*(cell[x,y].tree[t].location[2] - cell[x,y].tree[j].location[2]);
							end;
						end; //6
					end; //5
				end; //4
				
				if potmates_cnt = 0 then
				begin //4
					repeat //5
						ring := ring+1;

						xr := x+ring;
						yr := y+ring;
						
						l := 0;
						dir := 1;
						for i := 1 to ((ring*2-1)*4+4) do
						begin //6
							inc(l);
							
							xr := xr + directions[1,dir];
							yr := yr + directions[2,dir];
							
							hxr := maketorus(xr);
							hyr := maketorus(yr);
							
							if cell[hxr,hyr].no_trees > 0 then
							begin //7
								for m :=1 to cell[hxr,hyr].no_trees do
								begin //8
									if cell[hxr,hyr].tree[m].age > 1 then 
									begin //9
									
										inc(potmates_cnt);
										potmates[1,potmates_cnt] := hxr;
										potmates[2,potmates_cnt] := hyr;
										potmates[3,potmates_cnt] := m;	
										potmates_dist[potmates_cnt] := (cell[x,y].tree[t].location[1] - cell[hxr,hyr].tree[m].location[1])*(cell[x,y].tree[t].location[1] - cell[hxr,hyr].tree[m].location[1]) + (cell[x,y].tree[t].location[2] - cell[hxr,hyr].tree[m].location[2])*(cell[x,y].tree[t].location[2] - cell[hxr,hyr].tree[m].location[2]);
									
										stop := true;
									end; //9
								end; //8
							end; //7
							
							if l = ((ring*2-1)+1) then
							begin //7
							// direction change
							l :=0;
							dir := dir +1;
							end; //7
							
						end; //6
						
					until (stop = true);//5
				end; //4
				
				// mating, application of maternal costs, reproduction and dispersal
				if potmates_cnt > 0 then
				begin //4
					// now choose one matingpartner randomly from list
					// mchoice := random(potmates_cnt)+1;
					// choose nearest potential mating partner in euclidean distance
					// set mindist to first value
					mindist := potmates_dist[1];
					nearestneighbour := 1;
					if potmates_cnt > 1 then
					begin //5
						for pm := 2 to potmates_cnt do
						begin //6
							if potmates_dist[pm] < mindist then 
							begin //7
								mindist := potmates_dist[pm];
								nearestneighbour := pm;
							end;//7
						end;//6
					end;//5
					
					// use nearest neighbour as mate
					xf := potmates[1,nearestneighbour];
					yf := potmates[2,nearestneighbour];
					tf := potmates[3,nearestneighbour];

					// calculate mean and effective fertility including possible trade-off
					tradeOff_costs := 0;
					if dispCostCheckBox = 'maternal' then tradeOff_costs := mort_kernel(cell[x,y].tree[t].kernel);
					
					// determine effective no. of offspring
					no_offspring := poisson(fertility*(1-tradeOff_costs));
					
					//loop over all offspring
					if no_offspring > 0 then
					begin //5
						for o := 1 to no_offspring do
						begin //6
							// the new offspring disperses and eventually settles
							dispersal(x,y,t,xf,yf,tf);
						end; //6
					end;//5

				end //4
				else writeln('ERROR: no mates found!');
			end; //3
		end; //2
	end;
end;

//----------------------------------------------------------------------
// MAIN PROCEDURE ------------------------------------------------------
begin
	for run_no := 1 to RUN_MAX do
	begin
		// initialization
		initialize;

		for time:= 1 to TIME_MAX do
		begin // time loop
			write(run_no:10,'    ',time:10,'   ');
			// potential tree death: effects of age and competition
			death;
			//write(population_size:10,'    ');
			write(population_size_afterComp:10,'    ');
			// initilize
			mean_dispDist := 0;
			mean_dispDist_cnt := 0;
			// reproduction of trees and dispersal
			reproduction;
			if mean_dispDist_cnt > 0 then mean_dispDist := mean_dispDist/mean_dispDist_cnt
			else mean_dispDist := 999;
			writeln(mean_dispDist:10:2);
			
		end; // time loop
		
		// data analysis and output only once at the end of a simulation run
		if (time = TIME_MAX) then analysis;
	end;
end.
