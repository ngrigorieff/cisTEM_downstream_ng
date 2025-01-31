/*
 * scattering_potential.cpp

 *
 *  Created on: Oct 3, 2019
 *      Author: himesb
 */
#include "../../core/core_headers.h"
#include "scattering_potential.h"

ScatteringPotential::ScatteringPotential()
{

	int a = 5;

}

ScatteringPotential::~ScatteringPotential()
{
	if ( is_allocated_pdb_ensemble )
	{
		delete [] pdb_ensemble;
	}
}

void ScatteringPotential::InitPdbEnsemble(	float wanted_pixel_size, float do3d, int minimum_padding_x_and_y, int minimum_thickness_z,
											int max_number_of_noise_particles,
											float wanted_noise_particle_radius_as_mutliple_of_particle_radius,
											float wanted_noise_particle_radius_randomizer_lower_bound_as_praction_of_particle_radius,
											float wanted_noise_particle_radius_randomizer_upper_bound_as_praction_of_particle_radius,
											float wanted_tilt_angle_to_emulate)
{

	// backwards compatible with tigress where everything is double (ints would make more sense here.)
	long access_type_read = 0;
	long records_per_line = 1;

	pdb_ensemble = new PDB[number_of_pdbs] ;
	is_allocated_pdb_ensemble = true;
	// Initialize each of the PDB objects, this reads in and centers each PDB, but does not make any copies (instances) of the trajectories.

	for (int iPDB = 0; iPDB < number_of_pdbs ; iPDB++)
	{

//		// Need to find a way to pass this COM in. For now it was hardcoded.
//		if (do3d > 0.5)
//		{
//			// We've read through PDB[0] and got the center of mass. Overwrite this PDB with the given center of mass. Not sure why I put this in, but related to making an offset for LSU/SSU. FIXME
//			double COM[3] = {307.57,297.47,319.69};
//
//			pdb_ensemble[0] = PDB(pdb_file_names[iPDB],access_type_read, wanted_pixel_size, records_per_line, minimum_padding_x_and_y, minimum_thickness_z, COM);
//			number_of_pdbs = 1;
//		}
//		else
		{
			pdb_ensemble[iPDB] = PDB(pdb_file_names[iPDB],access_type_read, wanted_pixel_size, records_per_line, minimum_padding_x_and_y, minimum_thickness_z,
									 max_number_of_noise_particles,
									 wanted_noise_particle_radius_as_mutliple_of_particle_radius,
									 wanted_noise_particle_radius_randomizer_lower_bound_as_praction_of_particle_radius,
									 wanted_noise_particle_radius_randomizer_upper_bound_as_praction_of_particle_radius,
									 wanted_tilt_angle_to_emulate);
		}





	}
}

long ScatteringPotential::ReturnTotalNumberOfNonWaterAtoms()
{
	long number_of_non_water_atoms = 0;
	// Get a count of the total non water atoms
	for (int iPDB = 0; iPDB < number_of_pdbs; iPDB++)
	{
		//this->number_of_non_water_atoms += (pdb_ensemble[iPDB].number_of_atoms * this->particle_copy_number[iPDB]);
		number_of_non_water_atoms += (pdb_ensemble[iPDB].number_of_real_and_noise_atoms * pdb_ensemble[iPDB].number_of_particles_initialized);
	//		wxPrintf("%ld %d\n",pdb_ensemble[iPDB].number_of_atoms , pdb_ensemble[iPDB].number_of_particles_initialized);
		// These sizes will need to be determined by the min and max dimensions of the base shifted ensemble and removed from user input TODO

	}

	return number_of_non_water_atoms;
}



