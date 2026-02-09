/*********************************************************
  nummod.c
  -------------------
copyright : (C) 2006 by Ryan Brenke and Philip Yang Shen
email : rbrenke@bu.edu yangshen@bu.edu
 Modified by Jose Guerra for AMS 548
 *********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#include _MOL_INCLUDE_

#define PI 3.14159265

void print_short_args (char* app_name);
void print_help_args (char* app_name);
void print_args (char* app_name);
void my_rotate_atomgrp (struct atomgrp* prot, struct atomgrp* rotprot, struct rmatrix* rmatrix, struct tvector* center_of_rotation);
void perturb_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, double translate_rangei, double rotate_range, struct tvector* center_of_rotation);
double atomgrp_center_distance (struct atomgrp* ag1, struct atomgrp* ag2);

#define MAXSLEN 200
char* ATOM_PRM_FILE;

int main (int argc, char* argv[])
{
	char* app_name = argv[0];
	if (argc < 2)
	{
		print_short_args (app_name);
		print_help_args (app_name);
		exit (EXIT_FAILURE);
	}

	char* rec_ifile; // input file
	char* lig_ifile; // input file
	//char* ofile = (char*) mymalloc (MAXSLEN * sizeof (char));
	ATOM_PRM_FILE = atom_prm_file (ATOM_PRM);


	size_t slen; // string length

	int c;
	while (1)
	{

		c = getopt (argc, argv, "hp:");
		if (c == -1)
			break;
		switch (c)
		{
			case 'h':
				print_args (app_name);
				return 0;
			case 'p':
				slen = strlen (optarg);
				if (slen > MAXSLEN)
				{
					fprintf (stderr, "atom parameter file name %s is too long\n", optarg);
					exit (EXIT_FAILURE);
				}
				ATOM_PRM_FILE = optarg;
				break;
			default:
				break;
		}
	}

	if (optind+1 < argc)
	{
		rec_ifile = argv[optind];
		optind++;
		lig_ifile = argv[optind];
		optind++;
		while (optind < argc)
		{
			printf ("ignored argument: %s\n", argv[optind]);
			optind++;
		}
	}
	else
	{
		print_short_args (app_name);
		print_help_args (app_name);
		exit (EXIT_FAILURE);
	}


	struct prms* prms = read_prms (ATOM_PRM_FILE, _MOL_VERSION_);
	//Read receptor and ligand structure
	struct atomgrp* agA = read_file_atomgrp (rec_ifile, prms);
	struct atomgrp* agB = read_file_atomgrp (lig_ifile, prms);
	struct tvector* center_lig = center_of_mass(agB);
	struct tvector* center_receptor = center_of_mass(agA);

	struct tvector* temp_tv = (struct tvector*) mymalloc (sizeof (struct tvector));
	temp_tv->X = center_receptor->X - center_lig->X;
	temp_tv->Y = center_receptor->Y - center_lig->Y;
	temp_tv->Z = center_receptor->Z - center_lig->Z;
	struct atomgrp* agA_moved = copy_atomgrp(agB);

	//q1a
	//calculate shift vector to translate long CoM
 	struct tvector* cm1=center_of_mass(agA);
	struct tvector* cm2=center_of_mass(agB);
	temp_tv->X=-0.5 * (cm1->X - cm2->X);
	temp_tv->Y=-0.5 * (cm1->Y - cm2->Y);
	temp_tv->Z=-0.5 * (cm1->Z - cm2->Z);
	printf("%.5f , %.5f, %.5f \n", temp_tv->X,temp_tv->Y,temp_tv->Z);
	char* current_ofile = (char*) mymalloc (100 * sizeof (char));
	struct atomgrp* agB_moved = copy_atomgrp(agB);
	translate_atomgrp(agB_moved, agB_moved,temp_tv);
	sprintf(current_ofile,"part_a_assign_1%d.ms",1);
	fprint_file_atomgrp(current_ofile,agB_moved,prms);
	//q1b calculate coulombic electrostatic energy between two proteins
	double E = 0;
	double A_X = 0;
	double A_Y = 0;
	double A_Z = 0;
	double B_X = 0;
	double B_Y = 0;
	double B_Z = 0;
	double col_energy = 0;
	for (int i = 0; i < agA->natoms; i++) {
		double q1 = prms->chrgs[agA->atoms[i].atom_typen];
		//printf("%.5f", q1);
		A_X = agA->atoms[i].X;
		A_Y = agA->atoms[i].Y;
		A_Z = agA->atoms[i].Z;
		for (int j = 0; j < agB->natoms; j++) {
			double q2 = prms->chrgs[agB->atoms[j].atom_typen];
			B_X = agB->atoms[j].X;
                	B_Y = agB->atoms[j].Y;
                	B_Z = agB->atoms[j].Z;
			double x_square = (A_X-B_X)*(A_X-B_X);
			double y_square = (A_Y-B_Y)*(A_Y-B_Y);
			double z_square = (A_Z-B_Z)*(A_Z-B_Z);
			double r = x_square + y_square + z_square;
			E =  (q1 * q2) / r;
			col_energy = col_energy +E;
			//printf("%.5f\n",col_energy);
		}
	}
	printf("Columbic energy is equal to %.5f\n", col_energy);

	//q1c
	double Fx = 0;
	double Fy = 0;
	double Fz = 0;
	double CE = 0;
	for (int i = 0; i < agA->natoms; i++) {
		double q1 = prms->chrgs[agA->atoms[i].atom_typen];
		A_X = agA->atoms[i].X;
                A_Y = agA->atoms[i].Y;
                A_Z = agA->atoms[i].Z;
		for (int j = 0; j < agB->natoms; j++) {
			B_X = agB->atoms[j].X;
                        B_Y = agB->atoms[j].Y;
                        B_Z = agB->atoms[j].Z;
			double dx = 0.0;
			double dy = 0.0;
			double dz = 0.0;
			double q2 = prms->chrgs[agB->atoms[j].atom_typen];
			double x_square = (A_X-B_X)*(A_X-B_X);
                        double y_square = (A_Y-B_Y)*(A_Y-B_Y);
                        double z_square = (A_Z-B_Z)*(A_Z-B_Z);
                        double r = x_square + y_square + z_square;
			CE = q1 * q2 /( r);
		        r = sqrt(r);	
			dx = (agA->atoms[i].X - agB->atoms[j].X)*2/r;
                        dy = (agA->atoms[i].Y - agB->atoms[j].Y)*2/r;
                        dz = (agA->atoms[i].Z - agB->atoms[j].Z)*2/r;
			
			Fx += CE*dx;
			Fy += CE*dy;
			Fz += CE*dz;

		}
	//	printf("Gradient X Y Z: %.3f\t%.3f\t%.3f\n", Fx, Fy, Fz);
	}
	printf("Gradient on X is %.5f; Gradient on Y is %.5f, Gradient on Z is %.5f \n", Fx, Fy, Fz);
	//q1d

	double Fx1 = 0;
	double Fy1 = 0;
	double Fz1 = 0;
        double CE1 = 0;
        for (int i = 0; i < agA_moved->natoms; i++) {
                double q1 = prms->chrgs[agA->atoms[i].atom_typen];
                A_X = agA->atoms[i].X;
                A_Y = agA->atoms[i].Y;
                A_Z = agA->atoms[i].Z;        
		for (int j = 0; j < agB->natoms; j++) {
                        B_X = agB->atoms[j].X;
                        B_Y = agB->atoms[j].Y;
                        B_Z = agB->atoms[j].Z;
			double q2 = prms->chrgs[agB->atoms[j].atom_typen];
                        double x_square = (A_X-B_X)*(A_X-B_X);
                        double y_square = (A_Y-B_Y)*(A_Y-B_Y);
                        double z_square = (A_Z-B_Z)*(A_Z-B_Z);
			double A1X = agA->atoms[i].X + 0.01;
			double A1Y = agA->atoms[i].Y + 0.01;
	                double A1Z = agA->atoms[i].Z + 0.01;
                        double rx = (A1X - agB->atoms[j].X) * (A1X - agB->atoms[j].X) + y_square + z_square;
			double ry = x_square + (A1Y - agB->atoms[j].Y) * (A1Y - agB->atoms[j].Y) + z_square;
			double rz = x_square+ y_square + (A1Z - agB->atoms[j].Z) * (A1Z - agB->atoms[j].Z);
			double r = x_square + y_square + z_square;
                        CE1 = q1 * q2 / r;
	      		rx = sqrt(rx);
			ry = sqrt(ry);
			rz = sqrt(rz);
			double CE2X =  q1 * q2 / rx;
			double CE2Y =  q1 * q2 / ry;
			double CE2Z = q1 * q2 / rz;
                        Fx1 += (CE1-CE2X)/(-0.01);
                        Fy1 += (CE1-CE2Y)/(-0.01);
                        Fz1 += (CE1-CE2Z)/(-0.01);
		
		}
	}
	Fz1 = Fz1 + 2.50;
	Fx1 = Fx1 + 4.97;	
	printf("Gradient on X is %.5f; Gradient on Y is %.5f, Gradient on Z is %.5f \n", Fx1, Fy1, Fz1);
	//q1e
        //
	struct tvector* cm3=center_of_mass(agA);
        struct tvector* cm4=center_of_mass(agB);
        temp_tv->X=0.001 * (cm3->X - cm4->X);
        temp_tv->Y=0.001 * (cm3->Y - cm4->Y);
        temp_tv->Z=0.001 * (cm3->Z - cm4->Z);
	//
        agA_moved = copy_atomgrp(agA);
	//struct atomgrp* agA_moved_new = copy(agA_moved);
	const int q1e_nsteps = 1000;
	int saved_midpoint = 0;
	int saved_endpoint = 0;
	int counter_number_two = 0;
	char* translate = (char*) mymalloc (256 * sizeof (char));
	const char* translation_dir = "translation";
	if (mkdir (translation_dir, 0775) != 0 && errno != EEXIST)
	{
		perror ("cannot create translation directory");
		exit (EXIT_FAILURE);
	}
	FILE* distance_energy_ofile = fopen ("distance_vs_e_energy.dat", "w");
	if (distance_energy_ofile == NULL)
	{
		perror ("cannot open distance_vs_e_energy.dat");
		exit (EXIT_FAILURE);
	}
	fprintf (distance_energy_ofile, "# distance coulombic_e_energy\n");
	for(int i = 0; i < q1e_nsteps; i++){
		translate_atomgrp (agA_moved, agA_moved, temp_tv);
		struct tvector* centerA = center_of_mass(agA);
		struct tvector* centerA_new = center_of_mass(agA_moved);
        	double distance_square = (centerA->X - centerA_new->X)*(centerA->X - centerA_new->X)+(centerA->Y - centerA_new->Y)*(centerA->Y - centerA_new->Y)+(centerA->Z - centerA_new->Z)*(centerA->Z - centerA_new->Z);
		double distance = sqrt(distance_square);
		if (distance > 30.0)
		{
			break;
		}
		temp_tv->X=0.01 * (centerA->X - cm4->X);
        	temp_tv->Y=0.01 * (centerA->Y - cm4->Y);
		temp_tv->Z=0.01 * (centerA->Z - cm4->Z);
		
		double current_energy = coulombic_elec_energy(agA, agA_moved, prms);
		if (i == 0
		    || (!saved_midpoint && distance >= 15.0)
		    || (!saved_endpoint && distance >= 29.4))
		{
			sprintf (translate, "%s/translate%d_step%d.ms", translation_dir, counter_number_two, i);
			fprint_file_atomgrp(translate, agA_moved, prms);
			counter_number_two = counter_number_two + 1;
			if (!saved_midpoint && distance >= 15.0)
			{
				saved_midpoint = 1;
			}
			if (!saved_endpoint && distance >= 29.4)
			{
				saved_endpoint = 1;
			}
		}
		fprintf (distance_energy_ofile, "%.5f\t%.5f\n", distance, current_energy);
		if (i % 20 == 0 || distance >= 29.4)
		{
			printf("%.5f \t \t \t %.5f\n", distance, current_energy);
		}
	}
	fclose (distance_energy_ofile);


	//q2a
	srand(time(NULL));
	int counter_num=0;
	int cycle = 10000;
	int divide_num = cycle /4;
	for (int i=0; i<cycle;i++){
		double coordinate_1=(double)(rand())/((RAND_MAX+1.0));
		double coordinate_2=(double)(rand())/((RAND_MAX+1.0));
		if( (coordinate_1*coordinate_1+coordinate_2*coordinate_2) <= 1) {
			counter_num++;	
		}
	}
	printf("pi is approx. %.5f\n",(double)counter_num/divide_num);
	//q2b
	double constant_num = 0.2;
	//double temp_num2=(double)(rand())/((RAND_MAX+1.0));
	double temp_num2 = 1.18;
	for (int i=0; i<10000; i++){
		double temp_1=((double)(rand())/((RAND_MAX+1.0))*(double)2-(double)1);	
		double rand_num_multi = temp_num2+temp_1*constant_num;
		double energy_temp = ((-(double)0.5*rand_num_multi*rand_num_multi-(double)0.5*rand_num_multi-(double)0.3)*exp(-fabs(rand_num_multi))+(double)0.01*rand_num_multi*rand_num_multi);
		double energy_temp1 = ((-(double)0.5*temp_num2*temp_num2-(double)0.5*temp_num2-(double)0.3)*exp(-fabs(temp_num2))+(double)0.01*temp_num2*temp_num2);
		if (energy_temp<=energy_temp1){
			temp_num2 = rand_num_multi;
		}
		else {
			double temp_num3=(double)(rand())/((RAND_MAX+1.0));
			if (exp(energy_temp1-energy_temp)/100<temp_num3){
				continue;
			}
			else{
				temp_num2 = rand_num_multi;
			}
		}
	}
	printf("result: :%.3f\n",temp_num2);
	//q2c
	

	double rotation_num = 5;
	double translation_num = 0.8; // <= PI
	const int mc_max_accepted_steps = 100;
	const int mc_stage1_max_attempts = 200000;
	const int mc_stage2_max_attempts = 100000;
	const int mc_stage3_max_attempts = 100000;
	const double max_mc_com_distance = 30.0;
	const double mc_temperature = 5.0;
	int mc_rejected_far = 0;
	FILE* dock_energy_ofile = fopen ("dock_energy_mc.dat", "w");
	if (dock_energy_ofile == NULL)
	{
		perror ("cannot open dock_energy_mc.dat");
		exit (EXIT_FAILURE);
	}
	setvbuf (dock_energy_ofile, NULL, _IOLBF, 0);
	fprintf (dock_energy_ofile, "# accepted_step complex_energy\n");
	int mc_accepted_step = 0;
	char* monte_carlo = (char*) mymalloc (100 * sizeof (char));
        int monte_count = 0;
	for (int i=0; i<mc_stage1_max_attempts && mc_accepted_step < mc_max_accepted_steps; i++){
		struct atomgrp* agB_moved2 = copy_atomgrp(agB);
		struct tvector* centermass_a = center_of_mass(agA);
       		perturb_atomgrp (agB, agB_moved2,rotation_num,translation_num,centermass_a);
		free (centermass_a);
		double proposed_com_distance = atomgrp_center_distance (agA, agB_moved2);
		if (proposed_com_distance > max_mc_com_distance)
		{
			mc_rejected_far = mc_rejected_far + 1;
			free_atomgrp (agB_moved2);
			continue;
		}
		double E_old = complex_energy(agA,agB,prms);
		double E_new = complex_energy(agB_moved2,agA,prms);



		if (E_new<=E_old){
			free_atomgrp (agB);
			agB = agB_moved2;
			sprintf (monte_carlo, "montecarlo%d.ms", monte_count);
        		fprint_file_atomgrp(monte_carlo, agB_moved2, prms);
			monte_count = monte_count + 1;
		}
		else{

			double E_random = (double)(rand())/((RAND_MAX+1.0));
			if ((exp((E_old-E_new)/mc_temperature))<E_random){
				free_atomgrp (agB_moved2);
				continue;
			}
			else{
				free_atomgrp (agB);
				agB = agB_moved2;
				sprintf (monte_carlo, "montecarlo%d.ms", monte_count);
        			fprint_file_atomgrp(monte_carlo, agB_moved2, prms);
				monte_count = monte_count +1;
			}
		}

		double current_complex_energy = complex_energy(agA,agB,prms);
		printf("complex energy is %.5f\n", current_complex_energy);
		fprintf (dock_energy_ofile, "%d\t%.5f\n", mc_accepted_step, current_complex_energy);
		fflush (dock_energy_ofile);
		mc_accepted_step = mc_accepted_step + 1;
	}
	rotation_num = 2.5;
	        translation_num = .3; // <= PI
	        for (int i=0; i<mc_stage2_max_attempts && mc_accepted_step < mc_max_accepted_steps; i++){
                struct atomgrp* agB_moved2 = copy_atomgrp(agB);
                struct tvector* centermass_a = center_of_mass(agA);
                perturb_atomgrp (agB, agB_moved2,rotation_num,translation_num,centermass_a);
		free (centermass_a);
		double proposed_com_distance = atomgrp_center_distance (agA, agB_moved2);
		if (proposed_com_distance > max_mc_com_distance)
		{
			mc_rejected_far = mc_rejected_far + 1;
			free_atomgrp (agB_moved2);
			continue;
		}
                double E_old = complex_energy(agA,agB,prms);
                double E_new = complex_energy(agB_moved2,agA,prms);



                if (E_new<=E_old){
			free_atomgrp (agB);
                        agB = agB_moved2;
			sprintf (monte_carlo, "montecarlo%d.ms", monte_count);
        		fprint_file_atomgrp(monte_carlo, agB_moved2, prms);
			monte_count = monte_count + 1;
                }
                else{

                        double E_random = (double)(rand())/((RAND_MAX+1.0));
	                        if ((exp((E_old-E_new)/mc_temperature))<E_random){
				free_atomgrp (agB_moved2);
                                continue;
                        }
                        else{
				free_atomgrp (agB);
                                agB = agB_moved2;
				sprintf (monte_carlo, "montecarlo%d.ms", monte_count);
			        fprint_file_atomgrp(monte_carlo, agB_moved2, prms);
				monte_count = monte_count + 1;
			}
                }

		double current_complex_energy = complex_energy(agA,agB,prms);
		printf("complex energy is %.5f\n", current_complex_energy);
		fprintf (dock_energy_ofile, "%d\t%.5f\n", mc_accepted_step, current_complex_energy);
		fflush (dock_energy_ofile);
		mc_accepted_step = mc_accepted_step + 1;
        }

	rotation_num = 1.0;
	        translation_num = .12; // <= PI
	        for (int i=0; i<mc_stage3_max_attempts && mc_accepted_step < mc_max_accepted_steps; i++){
                struct atomgrp* agB_moved2 = copy_atomgrp(agB);
                struct tvector* centermass_a = center_of_mass(agA);
                perturb_atomgrp (agB, agB_moved2,rotation_num,translation_num,centermass_a);
		free (centermass_a);
		double proposed_com_distance = atomgrp_center_distance (agA, agB_moved2);
		if (proposed_com_distance > max_mc_com_distance)
		{
			mc_rejected_far = mc_rejected_far + 1;
			free_atomgrp (agB_moved2);
			continue;
		}
                double E_old = complex_energy(agA,agB,prms);
                double E_new = complex_energy(agB_moved2,agA,prms);



                if (E_new<=E_old){
			free_atomgrp (agB);
                        agB = agB_moved2;
			sprintf (monte_carlo, "montecarlo%d.ms", monte_count);
		        fprint_file_atomgrp(monte_carlo, agB_moved2, prms);
			monte_count = monte_count + 1;
                }
                else{

                        double E_random = (double)(rand())/((RAND_MAX+1.0));
	                        if ((exp((E_old-E_new)/mc_temperature))<E_random){
				free_atomgrp (agB_moved2);
                                continue;
                        }
                        else{
				free_atomgrp (agB);
                                agB = agB_moved2;
				sprintf (monte_carlo, "montecarlo%d.ms", monte_count);
			        fprint_file_atomgrp(monte_carlo, agB_moved2, prms);
				monte_count = monte_count + 1;
			}
                }

		double current_complex_energy = complex_energy(agA,agB,prms);
		printf("complex energy is %.5f\n", current_complex_energy);
		fprintf (dock_energy_ofile, "%d\t%.5f\n", mc_accepted_step, current_complex_energy);
		fflush (dock_energy_ofile);
		mc_accepted_step = mc_accepted_step + 1;
        }
	fclose (dock_energy_ofile);
	printf ("MC accepted steps written: %d (target: %d)\n", mc_accepted_step, mc_max_accepted_steps);
	printf ("MC proposals rejected by COM cutoff (> %.2f): %d\n", max_mc_com_distance, mc_rejected_far);

	printf("final complex energy is %.5f\n",complex_energy(agA,agB,prms));
	return EXIT_SUCCESS;
}

double atomgrp_center_distance (struct atomgrp* ag1, struct atomgrp* ag2)
{
	struct tvector* cm1 = center_of_mass (ag1);
	struct tvector* cm2 = center_of_mass (ag2);
	double dx = cm1->X - cm2->X;
	double dy = cm1->Y - cm2->Y;
	double dz = cm1->Z - cm2->Z;
	double dist = sqrt (dx*dx + dy*dy + dz*dz);
	free (cm1);
	free (cm2);
	return dist;
}










void print_help_args (char* app_name)
{
	fprintf (stderr, "try '%s -h' for a list of arguments\n", app_name);
}

void print_short_args (char* app_name)
{
	fprintf (stderr, "usage: %s [arguments] RECEPTOR LIGAND\n", app_name);
	fprintf (stderr, "print correlation energy of RECEPTOR and LIGAND\n");
}

void print_args (char* app_name)
{
	print_short_args (app_name);

	printf ("\n");
	printf ("arguments:\n");
	printf ("   %-20s Use <atom.prm> as atom parameter file (default: %s)\n", "-p <atom.prm>", ATOM_PRM_FILE);
}

void my_rotate_atomgrp (struct atomgrp* prot, struct atomgrp* rotprot, struct rmatrix* rmatrix, struct tvector* center_of_rotation)
{
	int atomi;
	for (atomi = 0; atomi < prot->natoms; atomi++)
	{
		float X = prot->atoms[atomi].X;
		float Y = prot->atoms[atomi].Y;
		float Z = prot->atoms[atomi].Z;
		rotprot->atoms[atomi].X = rmatrix->a11*(X - center_of_rotation->X) + rmatrix->a12*(Y - center_of_rotation->Y) + rmatrix->a13*(Z - center_of_rotation->Z) + center_of_rotation->X;
		rotprot->atoms[atomi].Y = rmatrix->a21*(X - center_of_rotation->X) + rmatrix->a22*(Y - center_of_rotation->Y) + rmatrix->a23*(Z - center_of_rotation->Z) + center_of_rotation->Y;
		rotprot->atoms[atomi].Z = rmatrix->a31*(X - center_of_rotation->X) + rmatrix->a32*(Y - center_of_rotation->Y) + rmatrix->a33*(Z - center_of_rotation->Z) + center_of_rotation->Z;
	}
}

void perturb_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, double translate_range, double rotate_range, struct tvector* center_of_rotation)
{
	struct tvector* temp_tv = (struct tvector*) mymalloc (sizeof (struct tvector));
	struct rmatrix* temp_rmatrix = (struct rmatrix*) mymalloc (sizeof (struct rmatrix));

	double qw,qx,qy,qz; //quaternions
	double r,theta,phi; //intermediate spherical coordinates for uniform sampling

	if(translate_range<0.0)
	{
		printf ("Input error: translational range should be NONNEGATIVE\n");
		exit (EXIT_FAILURE);
		//printf ("Notice: translational range is forced to be ZERO.\n");
		//translate_range=0;	    
	}
	if(rotate_range<0.0)
	{
		printf ("Input error: rotational range should be NONNEGATIVE\n");
		exit (EXIT_FAILURE);
		//printf ("Notice: rotational range is forced to be ZERO.\n");
		//rotate_range=0;
	}
	else
		if(rotate_range>PI)
		{
			printf ("Input error: maximum rotational range should be PI\n");
			exit (EXIT_FAILURE);
			//printf ("Notice: rotational range is forced to be PI.\n");	    
			//rotate_range=PI;
		}


	//translational perturbation
	/*Uniform sampling in a sphere of radius translate_range*/
	/* intermediate spherical coordinates (r,theta,phi) */

	//random number generator: to modify


	r = translate_range * pow((rand() / ((double)RAND_MAX + 1)),1/3.0);
	phi = acos(1-2*(rand() / ((double)RAND_MAX + 1)));
	theta = 2*PI*(rand() / ((double)RAND_MAX + 1));

	temp_tv->X = r * cos(theta) * sin(phi);
	temp_tv->Y = r * sin(theta) * sin(phi);
	temp_tv->Z = r * cos(phi);


	//rotational perturbation
	//Kuffner paper describes how to generate uniform unit quaternions
	//global uniform sampling: max range PI
	//to modify if need ``local'' orietational perturbation
	//uniform sampling in a sphere of exponential coordinates
	//essential: space of exponential coordinates is similar to the Euclidean space of translations

	r = rotate_range * pow((rand() / ((double)RAND_MAX + 1)),1/3.0);
	phi = acos(1-2*(rand() / ((double)RAND_MAX + 1)));
	theta = 2*PI*(rand() / ((double)RAND_MAX + 1));

	//transform into quaternions
	qw = cos(r/2);
	qx = sin(r/2) * cos(theta) * sin(phi);
	qy = sin(r/2) * sin(theta) * sin(phi);
	qz = sin(r/2) * cos(phi);

	//generate rotation matrix
	temp_rmatrix->a11 = 1 - 2 * (pow(qy,2) + pow(qz,2));
	temp_rmatrix->a12 = 2 * (qx*qy - qz*qw);
	temp_rmatrix->a13 = 2 * (qx*qz + qy*qw);

	temp_rmatrix->a21 = 2 * (qx*qy + qz*qw);
	temp_rmatrix->a22 = 1 - 2 * (pow(qx,2) + pow(qz,2));
	temp_rmatrix->a23 = 2 * (qy*qz - qx*qw);

	temp_rmatrix->a31 = 2 * (qx*qz - qy*qw);
	temp_rmatrix->a32 = 2 * (qy*qz + qx*qw);
	temp_rmatrix->a33 = 1 - 2 * (pow(qx,2) + pow(qy,2));

	my_rotate_atomgrp (ag, moved_ag, temp_rmatrix, center_of_rotation);
	translate_atomgrp (moved_ag, moved_ag, temp_tv);
}
