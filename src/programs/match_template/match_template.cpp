#include "../../core/core_headers.h"


// Values for data that are passed around in the results.
const int number_of_output_images = 8; //mip, scaledmip, psi, theta, phi, pixel, defocus, sums, sqsums
const int number_of_meta_data_values = 6; // img_x, img_y, number cccs, histogram values.
const int MAX_ALLOWED_NUMBER_OF_PEAKS = 1000; // An error will be thrown and job aborted if this number of peaks is exceeded in the make template results block
const int blur_hack_n_frames = 64;
class AggregatedTemplateResult
{
public:


	int image_number;
	int number_of_received_results;
	float total_number_of_ccs;

	float *collated_data_array;
	float *collated_mip_data;
	float *collated_psi_data;
	float *collated_theta_data;
	float *collated_phi_data;
	float *collated_defocus_data;
	float *collated_pixel_size_data;
	float *collated_pixel_sums;
	float *collated_pixel_square_sums;
	long *collated_histogram_data;

	AggregatedTemplateResult();
	~AggregatedTemplateResult();
	void AddResult(float *result_array, long array_size, int result_number, int number_of_expected_results);
};




WX_DECLARE_OBJARRAY(AggregatedTemplateResult, ArrayOfAggregatedTemplateResults);
#include <wx/arrimpl.cpp> // this is a magic incantation which must be done!
WX_DEFINE_OBJARRAY(ArrayOfAggregatedTemplateResults);



// nasty globals to track histogram size

int histogram_number_of_points = 512;
float histogram_min = -12.5f;
float histogram_max = 22.5f;

class
MatchTemplateApp : public MyApp
{
	public:

	bool DoCalculation();
	void DoInteractiveUserInput();
	void MasterHandleProgramDefinedResult(float *result_array, long array_size, int result_number, int number_of_expected_results);
	void ProgramSpecificInit();

	// for master collation

	ArrayOfAggregatedTemplateResults aggregated_results;
	bool 		is_rotated_by_90 = false;

	float GetMaxJobWaitTimeInSeconds() {return 120.0f;}

	private:
};

class ImageProjectionComparison
{
public:
	Particle					*particle;
	ReconstructedVolume			*reference_volume;
	Image						*projection_image;
//	Image						*temp_image;
};

// This is the function which will be minimized
float FrealignObjectiveFunction(void *scoring_parameters, float *array_of_values)
{
	ImageProjectionComparison *comparison_object = reinterpret_cast < ImageProjectionComparison *> (scoring_parameters);
	comparison_object->particle->temp_parameters = comparison_object->particle->current_parameters;
	comparison_object->particle->UnmapParameters(array_of_values);

//	comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image, *comparison_object->particle->ctf_image,
//			comparison_object->particle->alignment_parameters, 0.0, 0.0,
//			comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true);

	if (comparison_object->particle->no_ctf_weighting) comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
			*comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
			comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, false, false, false, false);
	// Case for normal parameter refinement with weighting applied to particle images and 3D reference
	else if (comparison_object->particle->includes_reference_ssnr_weighting) comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
			*comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
			comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true, true, false, false);
	// Case for normal parameter refinement with weighting applied only to particle images
	else comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
			*comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
			comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true, false, true, true);

//	if (comparison_object->particle->origin_micrograph < 0) comparison_object->particle->origin_micrograph = 0;
//	comparison_object->particle->origin_micrograph++;
//	for (int i = 0; i < comparison_object->projection_image->real_memory_allocated; i++) {comparison_object->projection_image->real_values[i] *= fabs(comparison_object->projection_image->real_values[i]);}
//	comparison_object->projection_image->ForwardFFT();
//	comparison_object->projection_image->CalculateCrossCorrelationImageWith(comparison_object->particle->particle_image);
//	comparison_object->projection_image->SwapRealSpaceQuadrants();
//	comparison_object->projection_image->BackwardFFT();
//	comparison_object->projection_image->QuickAndDirtyWriteSlice("proj.mrc", comparison_object->particle->origin_micrograph);
//	comparison_object->projection_image->SwapRealSpaceQuadrants();
//	comparison_object->particle->particle_image->SwapRealSpaceQuadrants();
//	comparison_object->particle->particle_image->BackwardFFT();
//	comparison_object->particle->particle_image->QuickAndDirtyWriteSlice("part.mrc", comparison_object->particle->origin_micrograph);
//	comparison_object->particle->particle_image->SwapRealSpaceQuadrants();
//	exit(0);

//	float score =  	- comparison_object->particle->particle_image->GetWeightedCorrelationWithImage(*comparison_object->projection_image, comparison_object->particle->bin_index,
//			  comparison_object->particle->pixel_size / comparison_object->particle->signed_CC_limit)
//			- comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_float);
//	wxPrintf("psi, theta, phi, x, y, = %g, %g, %g, %g, %g, score = %g\n",
//			comparison_object->particle->alignment_parameters.ReturnPsiAngle(),
//			comparison_object->particle->alignment_parameters.ReturnThetaAngle(),
//			comparison_object->particle->alignment_parameters.ReturnPhiAngle(),
//			comparison_object->particle->alignment_parameters.ReturnShiftX(),
//			comparison_object->particle->alignment_parameters.ReturnShiftY(), score);
//	return score;
//	wxPrintf("sigma_noise, mask_volume, penalty = %g %g %g\n", comparison_object->particle->sigma_noise, comparison_object->particle->mask_volume,
//			comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_float));
	return 	- comparison_object->particle->particle_image->GetWeightedCorrelationWithImage(*comparison_object->projection_image, comparison_object->particle->bin_index,
			  comparison_object->particle->pixel_size / comparison_object->particle->signed_CC_limit)
			- comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_parameters);
		// This penalty term assumes a Gaussian x,y distribution that is probably not correct in most cases. It might be better to leave it out.

}

IMPLEMENT_APP(MatchTemplateApp)

void MatchTemplateApp::ProgramSpecificInit()
{
}

// override the DoInteractiveUserInput

void MatchTemplateApp::DoInteractiveUserInput()
{
	wxString	input_search_images;
	wxString	input_reconstruction;

	wxString    mip_output_file;
	wxString    best_psi_output_file;
	wxString    best_theta_output_file;
	wxString    best_phi_output_file;
	wxString    best_defocus_output_file;
	wxString	best_pixel_size_output_file;

	wxString    output_histogram_file;
	wxString    correlation_std_output_file;
	wxString    correlation_avg_output_file;
	wxString    scaled_mip_output_file;

	float		pixel_size = 1.0f;
	float		voltage_kV = 300.0f;
	float		spherical_aberration_mm = 2.7f;
	float		amplitude_contrast = 0.07f;
	float 		defocus1 = 10000.0f;
	float		defocus2 = 10000.0f;;
	float		defocus_angle;
	float 		phase_shift;
	float		low_resolution_limit = 300.0;
	float		high_resolution_limit = 8.0;
	float		angular_step = 5.0;
	int			best_parameters_to_keep = 20;
	float 		defocus_search_range = 500;
	float 		defocus_step = 50;
	float 		pixel_size_search_range = 0.1f;
	float 		pixel_size_step = 0.02f;
	float		padding = 1.0;
	bool		ctf_refinement = false;
	float		particle_radius_angstroms = 0.0f;
	wxString	my_symmetry = "C1";
	float 		in_plane_angular_step = 0;
	bool 		use_gpu_input = false;
	int			max_threads = 1; // Only used for the GPU code

	UserInput *my_input = new UserInput("MatchTemplate", 1.00);

	input_search_images = my_input->GetFilenameFromUser("Input images to be searched", "The input image stack, containing the images that should be searched", "image_stack.mrc", true);
	input_reconstruction = my_input->GetFilenameFromUser("Input template reconstruction", "The 3D reconstruction from which projections are calculated", "reconstruction.mrc", true);
	mip_output_file = my_input->GetFilenameFromUser("Output MIP file", "The file for saving the maximum intensity projection image", "mip.mrc", false);
	scaled_mip_output_file = my_input->GetFilenameFromUser("Output Scaled MIP file", "The file for saving the maximum intensity projection image divided by correlation variance", "mip_scaled.mrc", false);
	best_psi_output_file = my_input->GetFilenameFromUser("Output psi file", "The file for saving the best psi image", "psi.mrc", false);
	best_theta_output_file = my_input->GetFilenameFromUser("Output theta file", "The file for saving the best psi image", "theta.mrc", false);
	best_phi_output_file = my_input->GetFilenameFromUser("Output phi file", "The file for saving the best psi image", "phi.mrc", false);
	best_defocus_output_file = my_input->GetFilenameFromUser("Output defocus file", "The file for saving the best defocus image", "defocus.mrc", false);
	best_pixel_size_output_file = my_input->GetFilenameFromUser("Output pixel size file", "The file for saving the best pixel size image", "pixel_size.mrc", false);
	correlation_avg_output_file = my_input->GetFilenameFromUser("Correlation average value", "The file for saving the average value of all correlation images", "corr_average.mrc", false);
	correlation_std_output_file = my_input->GetFilenameFromUser("Correlation variance output file", "The file for saving the variance of all correlation images", "corr_variance.mrc", false);
	output_histogram_file = my_input->GetFilenameFromUser("Output histogram of correlation values", "histogram of all correlation values", "histogram.txt", false);
	pixel_size = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
	voltage_kV = my_input->GetFloatFromUser("Beam energy (keV)", "The energy of the electron beam used to image the sample in kilo electron volts", "300.0", 0.0);
	spherical_aberration_mm = my_input->GetFloatFromUser("Spherical aberration (mm)", "Spherical aberration of the objective lens in millimeters", "2.7");
	amplitude_contrast = my_input->GetFloatFromUser("Amplitude contrast", "Assumed amplitude contrast", "0.07", 0.0, 1.0);
	defocus1 = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
	defocus2 = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
	defocus_angle = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
	phase_shift = my_input->GetFloatFromUser("Phase Shift (degrees)", "Additional phase shift in degrees", "0.0");
//	low_resolution_limit = my_input->GetFloatFromUser("Low resolution limit (A)", "Low resolution limit of the data used for alignment in Angstroms", "300.0", 0.0);
	high_resolution_limit = my_input->GetFloatFromUser("High resolution limit (A)", "High resolution limit of the data used for alignment in Angstroms", "8.0", 0.0);
	angular_step = my_input->GetFloatFromUser("Out of plane angular step (0.0 = set automatically)", "Angular step size for global grid search", "0.0", 0.0);
	in_plane_angular_step = my_input->GetFloatFromUser("In plane angular step (0.0 = set automatically)", "Angular step size for in-plane rotations during the search", "0.0", 0.0);
//	best_parameters_to_keep = my_input->GetIntFromUser("Number of top hits to refine", "The number of best global search orientations to refine locally", "20", 1);
	defocus_search_range = my_input->GetFloatFromUser("Defocus search range (A)", "Search range (-value ... + value) around current defocus", "500.0", 0.0);
	defocus_step = my_input->GetFloatFromUser("Defocus step (A) (0.0 = no search)", "Step size used in the defocus search", "50.0", 0.0);
	pixel_size_search_range = my_input->GetFloatFromUser("Pixel size search range (A)", "Search range (-value ... + value) around current pixel size", "0.1", 0.0);
	pixel_size_step = my_input->GetFloatFromUser("Pixel size step (A) (0.0 = no search)", "Step size used in the pixel size search", "0.01", 0.0);
	padding = my_input->GetFloatFromUser("Padding factor", "Factor determining how much the input volume is padded to improve projections", "1.0", 1.0, 2.0);
//	ctf_refinement = my_input->GetYesNoFromUser("Refine defocus", "Should the particle defocus be refined?", "No");
	particle_radius_angstroms = my_input->GetFloatFromUser("Mask radius for global search (A) (0.0 = max)", "Radius of a circular mask to be applied to the input images during global search", "0.0", 0.0);
	my_symmetry = my_input->GetSymmetryFromUser("Template symmetry", "The symmetry of the template reconstruction", "C1");
#ifdef ENABLEGPU
	use_gpu_input = my_input->GetYesNoFromUser("Use GPU", "Offload expensive calcs to GPU","No");
	max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#endif


	int first_search_position = -1;
	int last_search_position = -1;
	int image_number_for_gui = 0;
	int number_of_jobs_per_image_in_gui = 0;
	float min_peak_radius = 10.0f;

	wxString directory_for_results = "/dev/null"; // shouldn't be used in interactive
	wxString result_filename = "/dev/null"; // shouldn't be used in interactive

	delete my_input;


	my_current_job.ManualSetArguments("ttffffffffffifffffbfftttttttttftiiiitttfbi",	input_search_images.ToUTF8().data(),
															input_reconstruction.ToUTF8().data(),
															pixel_size,
															voltage_kV,
															spherical_aberration_mm,
															amplitude_contrast,
															defocus1,
															defocus2,
															defocus_angle,
															low_resolution_limit,
															high_resolution_limit,
															angular_step,
															best_parameters_to_keep,
															defocus_search_range,
															defocus_step,
															pixel_size_search_range,
															pixel_size_step,
															padding,
															ctf_refinement,
															particle_radius_angstroms,
															phase_shift,
															mip_output_file.ToUTF8().data(),
															best_psi_output_file.ToUTF8().data(),
															best_theta_output_file.ToUTF8().data(),
															best_phi_output_file.ToUTF8().data(),
															best_defocus_output_file.ToUTF8().data(),
															best_pixel_size_output_file.ToUTF8().data(),
															scaled_mip_output_file.ToUTF8().data(),
															correlation_std_output_file.ToUTF8().data(),
															my_symmetry.ToUTF8().data(),
															in_plane_angular_step,
															output_histogram_file.ToUTF8().data(),
															first_search_position,
															last_search_position,
															image_number_for_gui,
															number_of_jobs_per_image_in_gui,
															correlation_avg_output_file.ToUTF8().data(),
															directory_for_results.ToUTF8().data(),
															result_filename.ToUTF8().data(),
															min_peak_radius,
															use_gpu_input,
															max_threads);
}

// override the do calculation method which will be what is actually run..

bool MatchTemplateApp::DoCalculation()
{
	/*
	wxString input_particle_images 				= my_current_job.arguments[0].ReturnStringArgument(); // global
	wxString input_parameter_file 				= my_current_job.arguments[1].ReturnStringArgument(); // not sure
	wxString input_reconstruction				= my_current_job.arguments[2].ReturnStringArgument(); // global
	wxString input_reconstruction_statistics 	= my_current_job.arguments[3].ReturnStringArgument(); // global
	bool	 use_statistics						= my_current_job.arguments[4].ReturnBoolArgument();   // global
	wxString ouput_matching_projections 		= my_current_job.arguments[5].ReturnStringArgument(); // ignore (always false)
	wxString ouput_parameter_file				= my_current_job.arguments[6].ReturnStringArgument(); // not sure par file
	wxString ouput_shift_file					= my_current_job.arguments[7].ReturnStringArgument(); // not sure output
	wxString my_symmetry						= my_current_job.arguments[8].ReturnStringArgument(); // global
	int		 first_particle						= my_current_job.arguments[9].ReturnIntegerArgument(); // local (effectively ignore)
	int		 last_particle						= my_current_job.arguments[10].ReturnIntegerArgument(); // local (effectively ignore)
	float	 percent_used						= my_current_job.arguments[11].ReturnFloatArgument();
	float 	 pixel_size							= my_current_job.arguments[12].ReturnFloatArgument(); // local
	float    voltage_kV							= my_current_job.arguments[13].ReturnFloatArgument(); // local
	float 	 spherical_aberration_mm			= my_current_job.arguments[14].ReturnFloatArgument(); // local
	float    amplitude_contrast					= my_current_job.arguments[15].ReturnFloatArgument(); // local
	float	 molecular_mass_kDa					= my_current_job.arguments[16].ReturnFloatArgument(); // global
	float    inner_mask_radius					= my_current_job.arguments[17].ReturnFloatArgument(); // global
	float    outer_mask_radius					= my_current_job.arguments[18].ReturnFloatArgument(); // global
	float    low_resolution_limit				= my_current_job.arguments[19].ReturnFloatArgument(); // global
	float    high_resolution_limit				= my_current_job.arguments[20].ReturnFloatArgument(); // global
	float	 signed_CC_limit					= my_current_job.arguments[21].ReturnFloatArgument(); // global
	float	 classification_resolution_limit	= my_current_job.arguments[22].ReturnFloatArgument(); // global
	float    particle_radius_angstroms					= my_current_job.arguments[23].ReturnFloatArgument(); // global
	float	 high_resolution_limit_search		= my_current_job.arguments[24].ReturnFloatArgument(); // global
	float	 angular_step						= my_current_job.arguments[25].ReturnFloatArgument(); // global
	int		 best_parameters_to_keep			= my_current_job.arguments[26].ReturnIntegerArgument(); // global
	float	 max_search_x						= my_current_job.arguments[27].ReturnFloatArgument(); // global
	float	 max_search_y						= my_current_job.arguments[28].ReturnFloatArgument(); // global
	refine_particle.mask_center_2d_x			= my_current_job.arguments[29].ReturnFloatArgument(); // global
	refine_particle.mask_center_2d_y			= my_current_job.arguments[30].ReturnFloatArgument(); // global
	refine_particle.mask_center_2d_z			= my_current_job.arguments[31].ReturnFloatArgument(); // global
	refine_particle.mask_radius_2d				= my_current_job.arguments[32].ReturnFloatArgument(); // global
	float	 defocus_search_range				= my_current_job.arguments[33].ReturnFloatArgument(); // global
	float	 defocus_step						= my_current_job.arguments[34].ReturnFloatArgument(); // global
	float	 padding							= my_current_job.arguments[35].ReturnFloatArgument(); // global
//	float	 filter_constant					= my_current_job.arguments[35].ReturnFloatArgument();
	bool	 global_search						= my_current_job.arguments[36].ReturnBoolArgument(); // global
	bool	 local_refinement					= my_current_job.arguments[37].ReturnBoolArgument(); // global
// Psi, Theta, Phi, ShiftX, ShiftY
	refine_particle.parameter_map[3]			= my_current_job.arguments[38].ReturnBoolArgument(); //global
	refine_particle.parameter_map[2]			= my_current_job.arguments[39].ReturnBoolArgument(); //global
	refine_particle.parameter_map[1]			= my_current_job.arguments[40].ReturnBoolArgument(); // global
	refine_particle.parameter_map[4]			= my_current_job.arguments[41].ReturnBoolArgument(); // global
	refine_particle.parameter_map[5]			= my_current_job.arguments[42].ReturnBoolArgument(); // global
	bool 	 calculate_matching_projections		= my_current_job.arguments[43].ReturnBoolArgument(); // global - but ignore
	refine_particle.apply_2D_masking			= my_current_job.arguments[44].ReturnBoolArgument(); // global
	bool	 ctf_refinement						= my_current_job.arguments[45].ReturnBoolArgument(); // global
	bool	 normalize_particles				= my_current_job.arguments[46].ReturnBoolArgument();
	bool	 invert_contrast					= my_current_job.arguments[47].ReturnBoolArgument(); // global - but ignore.
	bool	 exclude_blank_edges				= my_current_job.arguments[48].ReturnBoolArgument();
	bool	 normalize_input_3d					= my_current_job.arguments[49].ReturnBoolArgument();
	bool	 threshold_input_3d					= my_current_job.arguments[50].ReturnBoolArgument();
	bool	 local_global_refine				= my_current_job.arguments[51].ReturnBoolArgument();
	int		 current_class						= my_current_job.arguments[52].ReturnIntegerArgument(); // global - but ignore.
	bool	 ignore_input_angles				= my_current_job.arguments[53].ReturnBoolArgument(); // during global search, ignore the starting parameters (this helps reduce bias)
*/

	wxDateTime start_time = wxDateTime::Now();

	wxString	input_search_images_filename = my_current_job.arguments[0].ReturnStringArgument();
	wxString	input_reconstruction_filename = my_current_job.arguments[1].ReturnStringArgument();
	float		pixel_size = my_current_job.arguments[2].ReturnFloatArgument();
	float		voltage_kV = my_current_job.arguments[3].ReturnFloatArgument();
	float		spherical_aberration_mm = my_current_job.arguments[4].ReturnFloatArgument();
	float		amplitude_contrast = my_current_job.arguments[5].ReturnFloatArgument();
	float 		defocus1 = my_current_job.arguments[6].ReturnFloatArgument();
	float		defocus2 = my_current_job.arguments[7].ReturnFloatArgument();
	float		defocus_angle = my_current_job.arguments[8].ReturnFloatArgument();;
	float		low_resolution_limit = my_current_job.arguments[9].ReturnFloatArgument();
	float		high_resolution_limit_search = my_current_job.arguments[10].ReturnFloatArgument();
	float		angular_step = my_current_job.arguments[11].ReturnFloatArgument();
	int			best_parameters_to_keep = my_current_job.arguments[12].ReturnIntegerArgument();
	float 		defocus_search_range = my_current_job.arguments[13].ReturnFloatArgument();
	float 		defocus_step = my_current_job.arguments[14].ReturnFloatArgument();
	float 		pixel_size_search_range = my_current_job.arguments[15].ReturnFloatArgument();
	float 		pixel_size_step = my_current_job.arguments[16].ReturnFloatArgument();
	float		padding = my_current_job.arguments[17].ReturnFloatArgument();
	bool		ctf_refinement = my_current_job.arguments[18].ReturnBoolArgument();
	float		particle_radius_angstroms = my_current_job.arguments[19].ReturnFloatArgument();
	float 		phase_shift = my_current_job.arguments[20].ReturnFloatArgument();
	wxString    mip_output_file = my_current_job.arguments[21].ReturnStringArgument();
	wxString    best_psi_output_file = my_current_job.arguments[22].ReturnStringArgument();
	wxString    best_theta_output_file = my_current_job.arguments[23].ReturnStringArgument();
	wxString    best_phi_output_file = my_current_job.arguments[24].ReturnStringArgument();
	wxString    best_defocus_output_file = my_current_job.arguments[25].ReturnStringArgument();
	wxString    best_pixel_size_output_file = my_current_job.arguments[26].ReturnStringArgument();
	wxString    scaled_mip_output_file = my_current_job.arguments[27].ReturnStringArgument();
	wxString    correlation_avg_output_file = my_current_job.arguments[28].ReturnStringArgument();
	wxString 	my_symmetry = my_current_job.arguments[29].ReturnStringArgument();
	float		in_plane_angular_step = my_current_job.arguments[30].ReturnFloatArgument();
	wxString    output_histogram_file = my_current_job.arguments[31].ReturnStringArgument();
	int 		first_search_position = my_current_job.arguments[32].ReturnIntegerArgument();
	int 		last_search_position = my_current_job.arguments[33].ReturnIntegerArgument();
	int 		image_number_for_gui = my_current_job.arguments[34].ReturnIntegerArgument();
	int 		number_of_jobs_per_image_in_gui = my_current_job.arguments[35].ReturnIntegerArgument();
	wxString    correlation_std_output_file = my_current_job.arguments[36].ReturnStringArgument();
	wxString	directory_for_results = my_current_job.arguments[37].ReturnStringArgument();
	wxString	result_output_filename = my_current_job.arguments[38].ReturnStringArgument();
	float		min_peak_radius = my_current_job.arguments[39].ReturnFloatArgument();
	bool		use_gpu   = my_current_job.arguments[40].ReturnBoolArgument();
	int			max_threads =  my_current_job.arguments[41].ReturnIntegerArgument();

	if (is_running_locally == false) max_threads = number_of_threads_requested_on_command_line; // OVERRIDE FOR THE GUI, AS IT HAS TO BE SET ON THE COMMAND LINE...

	if (max_threads > 1)
	{
		MyAssertTrue(use_gpu, "Using more than one thread only works in the GPU implementation\n.");
	}
	else
	{
		if(use_gpu)
		{
			wxPrintf("Warning, you are only using one thread on the GPU. Suggested minimum is 2. Check compute saturation using nvidia-smi -l 1\n");
		}
	}


	/*wxPrintf("input image = %s\n", input_search_images_filename);
	wxPrintf("input reconstruction= %s\n", input_reconstruction_filename);
	wxPrintf("pixel size = %f\n", pixel_size);
	wxPrintf("voltage = %f\n", voltage_kV);
	wxPrintf("Cs = %f\n", spherical_aberration_mm);
	wxPrintf("amp contrast = %f\n", amplitude_contrast);
	wxPrintf("defocus1 = %f\n", defocus1);
	wxPrintf("defocus2 = %f\n", defocus2);
	wxPrintf("defocus_angle = %f\n", defocus_angle);
	wxPrintf("low res limit = %f\n", low_resolution_limit);
	wxPrintf("high res limit = %f\n", high_resolution_limit_search);
	wxPrintf("angular step = %f\n", angular_step);
	wxPrintf("best params to keep = %i\n", best_parameters_to_keep);
	wxPrintf("defocus search range = %f\n", defocus_search_range);
	wxPrintf("defocus step = %f\n", defocus_step);
	wxPrintf("padding = %f\n", padding);
	wxPrintf("ctf_refinement = %i\n", int(ctf_refinement));
	wxPrintf("mask search radius = %f\n", mask_radius_search);
	wxPrintf("phase shift = %f\n", phase_shift);
	wxPrintf("symmetry = %s\n", my_symmetry);
	wxPrintf("in plane step = %f\n", in_plane_angular_step);
	wxPrintf("first location = %i\n", first_search_position);
	wxPrintf("last location = %i\n", last_search_position);
	*/

	ParameterMap parameter_map; // needed for euler search init
	//for (int i = 0; i < 5; i++) {parameter_map[i] = true;}
	parameter_map.SetAllTrue();

	float outer_mask_radius;
	float current_psi;
	float psi_step;
	float psi_max;
	float psi_start;
	float histogram_step;

	float expected_threshold;
	float actual_number_of_ccs_calculated;

	double histogram_min_scaled; // scaled for the x*y scaling which is only applied at the end.
	double histogram_step_scaled;// scaled for the x*y scaling which is only applied at the end.

	long *histogram_data;

	int current_bin;

	float temp_float;
	float variance;
	double temp_double;
	double temp_double_array[5];
	float factor_score;

	int number_of_rotations;
	long total_correlation_positions;
	long current_correlation_position;
	long total_correlation_positions_per_thread;
	long pixel_counter;

	int current_search_position;
	int current_x;
	int current_y;

	int factorizable_x;
	int factorizable_y;
	int factor_result_pos;
	int factor_result_neg;

	int defocus_i;
	int size_i;

	int i;

	long original_input_image_x;
	long original_input_image_y;
	int remove_npix_from_edge = 0;
	double sqrt_input_pixels;

	EulerSearch	global_euler_search;
	AnglesAndShifts angles;

	ImageFile input_search_image_file;
	ImageFile input_reconstruction_file;

	Curve whitening_filter;
	Curve number_of_terms;

	input_search_image_file.OpenFile(input_search_images_filename.ToStdString(), false);
	input_reconstruction_file.OpenFile(input_reconstruction_filename.ToStdString(), false);

	bool do_shift_blur_hack = false;
	bool do_exposure_filter_hack = false;
	float* shift_hack_x;
	float* shift_hack_y;
	float* shift_hack_d;

	if (do_shift_blur_hack)
	{
		wxPrintf("Doing the shift blur hack\n");
		wxString tmp_string = input_search_images_filename;

		bool is_set = false;
		if (tmp_string.AfterLast('/').StartsWith("May06_12.48.59_20_1"))shift_hack_x = new float[blur_hack_n_frames]{0.179864,0.226035,0.221800,0.101352,-0.028329,-0.082596,-0.054399,-0.002695,0.053346,0.122100,0.190697,0.215507,0.225322,0.271974,0.290356,0.289413,0.314659,0.307740,0.315842,0.341318,0.332869,0.323045,0.315857,0.289356,0.283353,0.282720,0.294014,0.296343,0.289211,0.306148,0.309726,0.305844,0.323079,0.316819,0.293598,0.282043,0.260122,0.241315,0.248519,0.246461,0.228426,0.232809,0.257193,0.268968,0.248897,0.225143,0.219199,0.245235,0.249922,0.214677,0.192343,0.219558,0.237092,0.214007,0.175625,0.145734,0.149653,0.159645,0.162509,0.165749,0.202934,0.221146,0.203262,0.206492} ; shift_hack_y = new float[blur_hack_n_frames]{-0.598946,-0.618100,-0.596117,-0.617086,-0.672865,-0.625347,-0.538435,-0.430179,-0.334469,-0.265841,-0.194734,-0.131746,-0.059495,0.013713,0.037070,0.048825,0.076130,0.097069,0.117404,0.117561,0.101573,0.096835,0.093847,0.102045,0.111480,0.093396,0.085382,0.101272,0.109509,0.122307,0.107102,0.104257,0.136126,0.126461,0.108804,0.117392,0.103071,0.106639,0.115067,0.077577,0.087286,0.092084,0.073823,0.089815,0.085058,0.073099,0.088257,0.086413,0.083106,0.084351,0.080511,0.096522,0.104486,0.097786,0.090173,0.101890,0.120414,0.129854,0.125678,0.112223,0.105365,0.088581,0.074687,0.047516} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.45.05_28_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.076964,-1.155506,-1.145443,-0.984338,-0.839038,-0.718226,-0.624298,-0.528425,-0.449828,-0.422612,-0.383440,-0.311139,-0.274371,-0.271291,-0.250447,-0.205356,-0.161654,-0.123947,-0.118827,-0.117010,-0.093291,-0.080379,-0.072606,-0.062191,-0.065524,-0.061106,-0.049512,-0.036031,-0.008826,0.000610,-0.011463,-0.013589,-0.005986,-0.001644,-0.005619,-0.028267,-0.042511,-0.020780,-0.001382,-0.020484,-0.040623,-0.040835,-0.019477,0.007584,0.011379,0.014358,0.036602,0.068958,0.093039,0.081050,0.037517,-0.001317,-0.021377,-0.010435,-0.017181,-0.023696,-0.005392,0.014390,0.041199,0.036931,0.020179,0.027502,0.028056,-0.000602} ; shift_hack_y = new float[blur_hack_n_frames]{-1.467819,-1.517214,-1.509818,-1.429056,-1.370687,-1.277726,-1.175693,-1.080080,-0.991035,-0.904557,-0.822457,-0.751113,-0.695840,-0.643655,-0.589811,-0.561242,-0.522274,-0.480975,-0.470147,-0.447988,-0.423892,-0.404792,-0.369750,-0.361115,-0.370128,-0.344448,-0.323633,-0.328483,-0.331914,-0.337692,-0.335106,-0.309707,-0.291854,-0.314094,-0.314378,-0.290137,-0.288601,-0.285002,-0.282582,-0.281342,-0.260333,-0.243799,-0.251907,-0.256573,-0.241887,-0.230431,-0.234336,-0.220039,-0.208872,-0.199933,-0.172255,-0.177624,-0.196024,-0.200350,-0.201365,-0.189013,-0.185747,-0.190471,-0.184814,-0.171504,-0.132875,-0.125272,-0.139523,-0.125905} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.50.26_29_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.860774,-1.962393,-1.958552,-1.754390,-1.585977,-1.428430,-1.260191,-1.087192,-0.923631,-0.793156,-0.710755,-0.632715,-0.553784,-0.502921,-0.452333,-0.412102,-0.385443,-0.331616,-0.275966,-0.241083,-0.205945,-0.198026,-0.192671,-0.153718,-0.153161,-0.165344,-0.158132,-0.172187,-0.165460,-0.162236,-0.183776,-0.167540,-0.131628,-0.129621,-0.126603,-0.123831,-0.119153,-0.106729,-0.106861,-0.119296,-0.127083,-0.112719,-0.090237,-0.075550,-0.064035,-0.084219,-0.106067,-0.108499,-0.130085,-0.148391,-0.158799,-0.181190,-0.169487,-0.140168,-0.116247,-0.096757,-0.100750,-0.105727,-0.093690,-0.100167,-0.117723,-0.118864,-0.097247,-0.081332} ; shift_hack_y = new float[blur_hack_n_frames]{-2.362815,-2.405574,-2.397419,-2.369353,-2.318621,-2.192277,-2.030024,-1.825586,-1.603813,-1.397378,-1.225573,-1.087353,-0.961576,-0.832365,-0.741603,-0.676229,-0.612156,-0.550275,-0.491437,-0.445003,-0.399288,-0.359874,-0.325914,-0.291133,-0.286091,-0.275162,-0.247801,-0.245983,-0.232563,-0.215429,-0.208775,-0.182715,-0.179121,-0.173521,-0.161846,-0.169603,-0.166114,-0.160631,-0.162945,-0.155614,-0.152739,-0.134762,-0.121702,-0.119578,-0.113578,-0.107618,-0.095590,-0.089334,-0.084852,-0.080421,-0.095136,-0.099723,-0.097732,-0.102186,-0.111258,-0.130107,-0.115832,-0.092733,-0.087592,-0.091352,-0.080471,-0.065658,-0.045364,-0.036452} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.34.54_26_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.925367,-1.023360,-1.029377,-0.789848,-0.553345,-0.422511,-0.338979,-0.295183,-0.284493,-0.263596,-0.259635,-0.243867,-0.243310,-0.243398,-0.217928,-0.220929,-0.234173,-0.240186,-0.247569,-0.247493,-0.251521,-0.245186,-0.233937,-0.215851,-0.174106,-0.172916,-0.158652,-0.104544,-0.075578,-0.073348,-0.068141,-0.074826,-0.067497,-0.062865,-0.103411,-0.136533,-0.102523,-0.073498,-0.068008,-0.045348,-0.045347,-0.044355,-0.022644,-0.047190,-0.080025,-0.079503,-0.082949,-0.081670,-0.068697,-0.065058,-0.065765,-0.051288,-0.060814,-0.073237,-0.054142,-0.053357,-0.085401,-0.086186,-0.071696,-0.070841,-0.065233,-0.073974,-0.058881,-0.019729} ; shift_hack_y = new float[blur_hack_n_frames]{-0.944960,-1.047369,-1.048236,-0.821022,-0.632156,-0.496834,-0.390889,-0.333466,-0.288562,-0.238561,-0.199403,-0.185534,-0.160070,-0.129960,-0.094711,-0.064952,-0.063925,-0.064606,-0.060546,-0.049025,-0.039550,-0.028910,-0.027580,-0.030181,-0.024559,-0.000989,-0.004457,-0.011646,0.001177,-0.005297,0.002212,0.000771,-0.005606,0.005308,0.007187,-0.004148,-0.023846,-0.022663,0.013564,0.007856,-0.024387,-0.028157,-0.022439,-0.000499,0.001265,-0.018326,-0.010222,0.005469,-0.007994,-0.023206,-0.017260,-0.014476,-0.011258,-0.008233,0.009394,0.036660,0.030886,0.020468,0.011494,-0.015163,-0.020307,-0.032419,-0.057815,-0.030260} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_15.16.42_47_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.052099,-1.104188,-1.106793,-0.994128,-0.875416,-0.776361,-0.688781,-0.629467,-0.567376,-0.471045,-0.408957,-0.362317,-0.319491,-0.294355,-0.266817,-0.232852,-0.207430,-0.194119,-0.168796,-0.151539,-0.143000,-0.119391,-0.128435,-0.157505,-0.153080,-0.140019,-0.119477,-0.107872,-0.103830,-0.074583,-0.057543,-0.052851,-0.052967,-0.067161,-0.073717,-0.065454,-0.064100,-0.050578,-0.043938,-0.061586,-0.065787,-0.065468,-0.056630,-0.056490,-0.076281,-0.071640,-0.052799,-0.035343,-0.011483,-0.038468,-0.058225,-0.048609,-0.045636,-0.054137,-0.064528,-0.070168,-0.050251,-0.031754,-0.034524,-0.024595,-0.004387,0.008977,0.005297,0.017223} ; shift_hack_y = new float[blur_hack_n_frames]{-0.851525,-0.911112,-0.900969,-0.798930,-0.706884,-0.603950,-0.538946,-0.473073,-0.371599,-0.292910,-0.249475,-0.206522,-0.163688,-0.133852,-0.112234,-0.091886,-0.091737,-0.074293,-0.064798,-0.061245,-0.025394,-0.001852,0.004550,0.027976,0.021355,0.000978,0.005850,-0.001426,-0.011208,-0.023253,-0.042684,-0.035370,-0.014224,-0.029532,-0.044951,-0.027068,-0.015952,-0.011837,-0.006072,-0.013504,-0.002449,0.007042,-0.011589,-0.004979,-0.010410,-0.020558,0.005220,-0.002825,-0.005116,0.019146,0.015376,0.019449,0.019553,-0.008157,-0.012305,-0.015200,-0.026907,-0.021098,-0.024895,-0.034552,-0.018035,0.000363,0.001768,-0.016115} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_18.59.55_91_1"))shift_hack_x = new float[blur_hack_n_frames]{0.335093,0.360048,0.366411,0.284682,0.208394,0.190741,0.197390,0.220488,0.235825,0.252992,0.290978,0.311947,0.327797,0.332588,0.332479,0.340972,0.334571,0.309322,0.289002,0.288453,0.274272,0.262773,0.237156,0.222483,0.227828,0.217993,0.199316,0.194724,0.194123,0.206117,0.219078,0.229439,0.235757,0.217369,0.216893,0.212926,0.193421,0.176284,0.154940,0.153922,0.160943,0.155043,0.138120,0.133777,0.132472,0.112589,0.102600,0.098796,0.081896,0.080267,0.089023,0.096501,0.108804,0.114018,0.125042,0.134081,0.147436,0.155262,0.159842,0.175957,0.173353,0.155282,0.151481,0.157438} ; shift_hack_y = new float[blur_hack_n_frames]{-0.622593,-0.604605,-0.587701,-0.716842,-0.838798,-0.797828,-0.662467,-0.500800,-0.364136,-0.202665,-0.066392,0.009897,0.064545,0.114971,0.159751,0.204639,0.223119,0.232300,0.256529,0.285616,0.285728,0.275097,0.276256,0.254215,0.235761,0.237031,0.231301,0.228926,0.233882,0.228425,0.237159,0.231567,0.205703,0.203736,0.204244,0.181518,0.156446,0.157795,0.170916,0.172400,0.163582,0.150901,0.147325,0.157373,0.138352,0.110806,0.110999,0.098351,0.071645,0.065782,0.066933,0.069055,0.085876,0.072648,0.053027,0.065741,0.061438,0.041545,0.027120,0.003312,0.001500,0.009583,0.002737,-0.003149} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_22.57.31_135_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.726345,-0.817400,-0.827650,-0.592464,-0.371410,-0.219880,-0.137395,-0.123083,-0.035964,0.011742,0.019388,0.062800,0.053864,0.057609,0.107990,0.116823,0.119373,0.143540,0.122102,0.137352,0.168638,0.141041,0.111969,0.115204,0.118926,0.105922,0.075467,0.053733,0.075469,0.065481,0.022264,-0.015956,-0.008829,-0.010989,-0.087616,-0.155410,-0.152974,-0.137751,-0.170691,-0.195398,-0.201594,-0.198194,-0.164020,-0.195483,-0.223812,-0.194127,-0.204971,-0.198835,-0.197780,-0.217313,-0.185827,-0.188606,-0.214686,-0.214619,-0.234431,-0.210187,-0.193519,-0.216066,-0.218227,-0.206242,-0.190538,-0.186030,-0.184307,-0.164987} ; shift_hack_y = new float[blur_hack_n_frames]{-1.403664,-1.553470,-1.574262,-1.188093,-0.766339,-0.526527,-0.402431,-0.355127,-0.260675,-0.175913,-0.166417,-0.171900,-0.194220,-0.200146,-0.141691,-0.094322,-0.093506,-0.054652,-0.033318,0.021073,0.081527,0.033530,-0.013745,0.001920,-0.015277,-0.014317,-0.030747,-0.072245,-0.053742,-0.077068,-0.118957,-0.137913,-0.133881,-0.141186,-0.228846,-0.269325,-0.251938,-0.279488,-0.294687,-0.318848,-0.356699,-0.320020,-0.274602,-0.322445,-0.342677,-0.330358,-0.345015,-0.339360,-0.322686,-0.327461,-0.333186,-0.319665,-0.318442,-0.315701,-0.312963,-0.314525,-0.304139,-0.276511,-0.280074,-0.280981,-0.262831,-0.252463,-0.239221,-0.247579} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_19.29.44_96_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.667405,-0.775482,-0.767187,-0.534745,-0.319124,-0.179742,-0.096038,-0.037986,0.021073,0.038239,0.029811,0.058417,0.064018,0.058926,0.084903,0.074887,0.070877,0.100696,0.084311,0.092148,0.114372,0.113263,0.111034,0.121374,0.118709,0.107642,0.087566,0.083473,0.102889,0.114671,0.105530,0.097021,0.091601,0.075670,0.078291,0.081250,0.054157,0.049573,0.065429,0.055789,0.079262,0.104325,0.086619,0.082216,0.083113,0.079628,0.109268,0.123697,0.112803,0.129963,0.148129,0.118812,0.107034,0.123410,0.110817,0.107371,0.100155,0.077868,0.085259,0.079799,0.055377,0.049039,0.038371,0.037731} ; shift_hack_y = new float[blur_hack_n_frames]{-0.736203,-0.882245,-0.873523,-0.546969,-0.271734,-0.081210,0.044140,0.104290,0.145846,0.147092,0.138860,0.159283,0.139472,0.121938,0.125230,0.108713,0.128651,0.136330,0.098322,0.093155,0.067704,0.065948,0.110936,0.102088,0.077039,0.071664,0.069981,0.105834,0.107241,0.071335,0.062435,0.071007,0.080476,0.082828,0.064348,0.055093,0.061786,0.043826,0.057181,0.050714,0.025324,0.050422,0.056140,0.020492,0.006206,-0.002589,-0.005492,0.004007,-0.010521,-0.007491,0.026286,0.032410,0.019299,0.028522,0.012375,0.000704,-0.014164,-0.046304,-0.036151,-0.028316,-0.043522,-0.054737,-0.062970,-0.046805} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_19.35.09_97_1"))shift_hack_x = new float[blur_hack_n_frames]{0.282880,0.286952,0.276426,0.299817,0.304548,0.279127,0.265679,0.235783,0.221279,0.201801,0.191177,0.194052,0.171191,0.134810,0.135471,0.142683,0.113660,0.079392,0.051446,0.072315,0.104166,0.083731,0.058721,0.073964,0.095828,0.111622,0.101914,0.081077,0.094408,0.094867,0.098807,0.084475,0.058882,0.059642,0.060420,0.034698,0.023660,0.021991,0.020496,0.034993,0.028614,0.009485,0.010734,0.008664,-0.000379,0.019107,0.027999,0.024220,0.026893,0.036400,0.052985,0.041243,-0.000188,-0.008819,-0.019317,-0.020796,-0.014962,-0.038841,-0.022291,0.015506,0.015190,0.020320,0.030451,0.026929} ; shift_hack_y = new float[blur_hack_n_frames]{0.208640,0.191309,0.200263,0.224292,0.227324,0.251435,0.284501,0.281041,0.302737,0.352398,0.399852,0.425843,0.432229,0.435573,0.450858,0.461590,0.440708,0.419098,0.398222,0.393627,0.408980,0.407247,0.393234,0.397549,0.378226,0.375515,0.369032,0.358562,0.365740,0.352843,0.319813,0.313680,0.306376,0.290370,0.282833,0.262121,0.259940,0.261719,0.247985,0.257384,0.263186,0.243683,0.229540,0.225254,0.225500,0.230358,0.218540,0.208234,0.192323,0.181736,0.184931,0.161210,0.142964,0.121733,0.094905,0.112671,0.129855,0.112574,0.110837,0.099646,0.099073,0.111367,0.094975,0.077630} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_20.34.46_108_1"))shift_hack_x = new float[blur_hack_n_frames]{1.358260,1.405086,1.399419,1.321078,1.248741,1.143784,1.049547,0.968444,0.875247,0.772522,0.688884,0.608526,0.532621,0.482482,0.454740,0.444263,0.426051,0.387376,0.370328,0.365753,0.342677,0.326084,0.325327,0.310464,0.294362,0.284012,0.276256,0.268254,0.236279,0.223606,0.214865,0.198512,0.196120,0.187294,0.174897,0.179711,0.151602,0.137608,0.163196,0.161847,0.152250,0.144606,0.144252,0.156396,0.152925,0.141814,0.141462,0.137212,0.126662,0.117160,0.110128,0.112115,0.107432,0.097239,0.089843,0.092037,0.087521,0.079041,0.074193,0.067565,0.072065,0.074454,0.060612,0.061257} ; shift_hack_y = new float[blur_hack_n_frames]{1.305334,1.330049,1.318686,1.307926,1.322193,1.284659,1.186285,1.061071,0.966502,0.884427,0.817917,0.759776,0.686114,0.639974,0.599183,0.561831,0.546582,0.505101,0.455352,0.420061,0.404221,0.403577,0.393753,0.368948,0.346134,0.335299,0.326989,0.304153,0.285684,0.251745,0.215217,0.214789,0.214229,0.203594,0.195267,0.183520,0.192156,0.195911,0.182829,0.180400,0.166828,0.150569,0.130814,0.101641,0.104693,0.108233,0.095714,0.104957,0.105836,0.110473,0.131534,0.113119,0.091962,0.088979,0.073921,0.067917,0.074746,0.064527,0.063445,0.072785,0.067999,0.067009,0.070272,0.068337} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_22.46.52_133_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.172645,-1.342608,-1.331696,-0.991606,-0.699194,-0.425268,-0.245793,-0.093280,0.040335,0.098677,0.148318,0.183739,0.204032,0.221149,0.203013,0.188542,0.214307,0.216180,0.197999,0.167308,0.157028,0.183773,0.186421,0.172240,0.185256,0.186262,0.185001,0.169554,0.164032,0.172675,0.166181,0.141329,0.121173,0.130064,0.154644,0.159859,0.176188,0.200665,0.191580,0.197048,0.218833,0.217226,0.186028,0.150831,0.120449,0.134813,0.135453,0.091232,0.094753,0.117514,0.110890,0.137573,0.170117,0.159636,0.184383,0.181551,0.167094,0.186168,0.158630,0.133877,0.157904,0.156737,0.162048,0.146995} ; shift_hack_y = new float[blur_hack_n_frames]{-0.949452,-1.067334,-1.054373,-0.816981,-0.572542,-0.381067,-0.283475,-0.204021,-0.135737,-0.078275,-0.037338,-0.027804,-0.017854,0.011109,0.015388,0.009809,-0.002214,-0.010735,0.016368,0.012873,0.004390,0.033485,0.063150,0.082577,0.073974,0.053951,0.055837,0.064823,0.051916,0.022926,0.008230,0.018948,0.002249,-0.012682,-0.005349,0.027855,0.023254,-0.007900,-0.018789,-0.006244,0.013670,0.000750,-0.016484,0.009745,0.023559,0.009565,0.005405,-0.014161,-0.014829,-0.031526,-0.048297,-0.037131,-0.015471,-0.031923,-0.036650,0.001057,0.008398,-0.022805,-0.027794,-0.012919,-0.003792,0.018394,-0.014471,-0.038928} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_16.17.41_59_1"))shift_hack_x = new float[blur_hack_n_frames]{-2.329161,-2.510585,-2.479475,-2.181293,-1.950779,-1.668302,-1.394380,-1.139286,-0.916945,-0.737014,-0.618840,-0.552155,-0.468922,-0.369114,-0.306351,-0.265349,-0.245306,-0.215485,-0.142964,-0.109324,-0.125118,-0.119350,-0.101063,-0.060571,-0.041776,-0.069729,-0.086292,-0.077120,-0.065230,-0.071450,-0.079295,-0.070887,-0.061287,-0.037153,-0.015752,-0.009026,0.000665,-0.011395,-0.031049,-0.024330,-0.016100,-0.021838,-0.029966,-0.033762,-0.029569,-0.014376,-0.009310,-0.026765,-0.015595,-0.011320,-0.007365,0.031451,0.047299,0.057634,0.063173,0.047888,0.065781,0.067123,0.039644,0.027443,0.011756,0.021399,0.034955,0.018568} ; shift_hack_y = new float[blur_hack_n_frames]{-6.238672,-6.631771,-6.540968,-5.956625,-5.546288,-4.916944,-4.145427,-3.392611,-2.851956,-2.367688,-1.975830,-1.654014,-1.359872,-1.160797,-0.979304,-0.829865,-0.695774,-0.555724,-0.454846,-0.373588,-0.307991,-0.249691,-0.193358,-0.175031,-0.158955,-0.129011,-0.100249,-0.055491,-0.044981,-0.025633,0.007832,0.005220,-0.007799,-0.002258,0.002514,0.001092,0.005298,0.022701,0.056633,0.089432,0.085794,0.069040,0.099860,0.117727,0.100071,0.094023,0.093721,0.100852,0.111989,0.104802,0.101039,0.095645,0.093565,0.106812,0.104707,0.093559,0.098871,0.102837,0.103501,0.115559,0.112015,0.110324,0.117097,0.119209} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_20.08.22_103_1"))shift_hack_x = new float[blur_hack_n_frames]{-2.653414,-3.371696,-3.016379,-1.983428,-1.923842,-1.694024,-1.426651,-1.077187,-0.846899,-2.109006,-3.257958,-3.044918,-2.926828,-2.836345,-1.431031,-1.361190,-1.251878,0.145692,0.116013,0.130698,0.194653,0.196211,0.186323,0.185209,0.174496,0.228351,0.254088,0.220209,0.239491,0.286109,0.335207,0.341416,0.332274,0.348861,0.369384,0.367614,0.344879,0.328191,0.361692,0.384744,0.378542,0.376136,0.376005,0.407344,0.399503,0.357515,0.335084,0.324961,0.332343,0.360090,0.363758,0.360550,0.374306,0.398763,0.412594,0.409267,0.391927,0.366823,0.365454,0.359523,0.340482,0.335603,0.328728,0.309971} ; shift_hack_y = new float[blur_hack_n_frames]{-2.558412,-3.163280,-2.851393,-2.034921,-2.047948,-1.840265,-1.577352,-1.207987,-0.936356,-1.921108,-2.831281,-2.652021,-2.524229,-2.373098,-1.154101,-1.140204,-1.086831,0.127418,0.159883,0.187764,0.229430,0.280537,0.319299,0.298260,0.272341,0.295991,0.288971,0.226572,0.212735,0.265477,0.298785,0.287354,0.292860,0.299984,0.316920,0.324155,0.289222,0.276864,0.294174,0.285897,0.282111,0.302859,0.315564,0.330343,0.333018,0.316426,0.317015,0.321163,0.307929,0.293351,0.268566,0.283529,0.317027,0.301408,0.290967,0.277563,0.288391,0.321789,0.298773,0.261520,0.261564,0.263696,0.261083,0.247183} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_20.13.40_104_1"))shift_hack_x = new float[blur_hack_n_frames]{0.083284,0.106872,0.111237,0.049423,-0.068108,-0.142842,-0.116819,-0.084634,-0.074175,-0.062029,-0.046381,-0.000347,0.057248,0.093066,0.127520,0.147441,0.143283,0.124464,0.111705,0.125548,0.111131,0.070572,0.064314,0.076053,0.084173,0.093399,0.080705,0.070077,0.064149,0.054178,0.054565,0.053553,0.033680,0.011666,0.010627,0.019668,0.009195,-0.010324,-0.008355,0.002646,0.015045,0.013636,0.013635,0.021111,0.023005,0.016877,0.016635,0.024547,0.015302,0.005181,0.021498,0.052316,0.068050,0.063480,0.038801,0.041668,0.043829,0.026880,0.028826,0.029800,0.035249,0.028112,0.031374,0.061257} ; shift_hack_y = new float[blur_hack_n_frames]{-1.009315,-0.998596,-0.984757,-1.086564,-1.237015,-1.239719,-1.130506,-0.973290,-0.799199,-0.664130,-0.571686,-0.466022,-0.355565,-0.251804,-0.179693,-0.164102,-0.123745,-0.109037,-0.104544,-0.074674,-0.074214,-0.086799,-0.075023,-0.059937,-0.048447,-0.036831,-0.034689,-0.040636,-0.040414,-0.036760,-0.033647,-0.005567,0.004663,-0.027079,-0.049057,-0.035679,-0.030826,-0.044848,-0.060938,-0.081367,-0.075155,-0.043607,-0.047732,-0.061207,-0.062581,-0.071299,-0.069798,-0.061498,-0.089949,-0.104430,-0.079430,-0.095956,-0.125889,-0.114057,-0.103747,-0.090962,-0.096492,-0.124774,-0.103709,-0.079170,-0.099224,-0.105729,-0.107438,-0.113966} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.55.43_30_1"))shift_hack_x = new float[blur_hack_n_frames]{-2.128584,-2.259270,-2.255443,-1.996976,-1.766649,-1.524842,-1.331799,-1.188035,-1.052832,-0.912493,-0.772897,-0.690421,-0.649310,-0.579585,-0.501470,-0.469086,-0.424069,-0.352555,-0.318817,-0.314178,-0.288765,-0.261919,-0.228214,-0.227436,-0.229851,-0.210923,-0.164262,-0.142239,-0.151886,-0.130606,-0.131153,-0.135931,-0.115138,-0.114740,-0.110843,-0.076165,-0.052329,-0.035038,-0.036957,-0.041238,-0.032478,-0.033194,-0.038658,-0.068656,-0.061372,-0.023526,-0.061001,-0.076647,-0.035624,-0.046560,-0.050116,-0.051245,-0.080593,-0.064247,-0.046556,-0.078784,-0.072332,-0.069879,-0.105316,-0.083305,-0.066170,-0.097722,-0.097110,-0.084330} ; shift_hack_y = new float[blur_hack_n_frames]{-2.394580,-2.480039,-2.445298,-2.343915,-2.325697,-2.185215,-1.959916,-1.764920,-1.604015,-1.425106,-1.240183,-1.071856,-0.910645,-0.806177,-0.728031,-0.635456,-0.605760,-0.571895,-0.516694,-0.502584,-0.448263,-0.365870,-0.328982,-0.275241,-0.226550,-0.229846,-0.217000,-0.220454,-0.236683,-0.221587,-0.206162,-0.210654,-0.181296,-0.165381,-0.173041,-0.164377,-0.162006,-0.170689,-0.174669,-0.170386,-0.140034,-0.108063,-0.120103,-0.116569,-0.112690,-0.104390,-0.084379,-0.084495,-0.099385,-0.080932,-0.062085,-0.043521,-0.029914,-0.042576,-0.041698,-0.017959,0.007625,0.022868,0.026691,0.023826,0.023458,0.019327,0.005548,-0.000084} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_22.23.11_130_1"))shift_hack_x = new float[blur_hack_n_frames]{0.531525,0.580900,0.570100,0.474237,0.445523,0.387184,0.324620,0.313295,0.270564,0.231598,0.230292,0.225128,0.223274,0.236373,0.212476,0.183915,0.186425,0.188811,0.159008,0.148563,0.122378,0.076512,0.081944,0.121732,0.119086,0.117967,0.109141,0.075757,0.095366,0.106095,0.076052,0.053818,0.042713,0.051935,0.081353,0.076095,0.074139,0.068617,0.067611,0.085233,0.079234,0.057625,0.046920,0.046103,0.049166,0.039029,0.016342,0.032132,0.058850,0.057345,0.048529,0.046060,0.083045,0.087214,0.039946,0.029266,0.033621,0.049004,0.055615,0.027255,0.020626,0.037953,0.046403,0.068638} ; shift_hack_y = new float[blur_hack_n_frames]{-0.700358,-0.670385,-0.658037,-0.781746,-0.897382,-0.889902,-0.823684,-0.767489,-0.674294,-0.557869,-0.485707,-0.428138,-0.385565,-0.340036,-0.298965,-0.270593,-0.257583,-0.248877,-0.223533,-0.181667,-0.174255,-0.180557,-0.188060,-0.168631,-0.140680,-0.140087,-0.149793,-0.141409,-0.111016,-0.097309,-0.097467,-0.094209,-0.100748,-0.101159,-0.110603,-0.123103,-0.122699,-0.119229,-0.094610,-0.089538,-0.096267,-0.095135,-0.094497,-0.070466,-0.060076,-0.098434,-0.095228,-0.055778,-0.046677,-0.056412,-0.040469,-0.050671,-0.041898,-0.039217,-0.085993,-0.051416,-0.017254,-0.062852,-0.088562,-0.091880,-0.076579,-0.070870,-0.116616,-0.106046} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_23.12.42_138_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.894239,-0.881506,-0.864724,-0.975693,-1.113738,-1.122720,-1.023114,-0.891643,-0.772309,-0.638840,-0.524054,-0.435592,-0.345563,-0.256818,-0.195706,-0.147018,-0.105098,-0.071504,-0.073259,-0.056623,-0.028229,-0.013456,-0.012503,-0.033131,-0.026938,0.023415,0.037512,0.021495,0.026402,0.031806,0.049032,0.054478,0.053478,0.052308,0.036279,0.030439,0.049182,0.048124,0.036679,0.016890,0.023143,0.067494,0.081567,0.073573,0.068669,0.061958,0.061542,0.054122,0.031997,0.002851,0.002165,0.027087,0.033516,0.028024,0.036813,0.054722,0.075620,0.069744,0.059008,0.057311,0.048819,0.051432,0.062210,0.065769} ; shift_hack_y = new float[blur_hack_n_frames]{-1.952593,-1.965846,-1.925356,-2.040171,-2.228294,-2.176815,-1.914235,-1.622159,-1.420263,-1.214512,-0.982982,-0.771413,-0.631448,-0.543916,-0.427556,-0.294835,-0.209375,-0.152004,-0.140127,-0.084797,-0.019460,-0.024034,-0.030591,-0.021870,-0.019026,-0.010321,-0.012340,-0.006628,0.017169,0.018943,0.018284,0.010983,0.034833,0.054872,0.046460,0.042243,0.028472,0.038640,0.060694,0.038329,0.038066,0.042139,0.031631,0.047476,0.040375,0.036527,0.046389,0.043559,0.071109,0.079764,0.064145,0.064204,0.052051,0.046578,0.040518,0.014917,0.006725,0.026656,0.043059,0.040634,0.055956,0.078972,0.074034,0.076813} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_16.58.30_67_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.170469,-1.210852,-1.196806,-1.154858,-1.132164,-1.032145,-0.897938,-0.810500,-0.729925,-0.609118,-0.518888,-0.435432,-0.375085,-0.347178,-0.291034,-0.239878,-0.206070,-0.170992,-0.162226,-0.144863,-0.133453,-0.130393,-0.102933,-0.111189,-0.119491,-0.104522,-0.107538,-0.082150,-0.060697,-0.082705,-0.096241,-0.085714,-0.072813,-0.065765,-0.077878,-0.107223,-0.104839,-0.082162,-0.089981,-0.117752,-0.124654,-0.101426,-0.072929,-0.074042,-0.080481,-0.058344,-0.037147,-0.029157,-0.029059,-0.028632,-0.019034,-0.014289,-0.017748,-0.006375,0.004973,-0.007154,-0.009563,0.005236,0.012431,0.010105,-0.003686,-0.003324,0.008172,0.016353} ; shift_hack_y = new float[blur_hack_n_frames]{-2.724787,-2.766529,-2.736019,-2.758398,-2.781773,-2.616501,-2.420876,-2.193379,-1.921038,-1.671321,-1.475414,-1.298691,-1.148922,-1.028029,-0.910824,-0.816648,-0.773145,-0.730381,-0.664624,-0.606776,-0.549402,-0.514371,-0.522937,-0.492499,-0.452197,-0.464844,-0.472497,-0.450839,-0.444025,-0.444502,-0.448862,-0.458959,-0.434231,-0.423699,-0.426918,-0.419390,-0.406562,-0.380095,-0.347098,-0.337752,-0.338418,-0.340582,-0.335156,-0.306736,-0.283265,-0.282046,-0.290739,-0.274339,-0.283855,-0.315079,-0.315734,-0.313916,-0.320155,-0.297769,-0.297430,-0.298881,-0.266894,-0.256896,-0.262123,-0.254775,-0.257391,-0.256016,-0.229822,-0.213800} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000} ; is_set = true;









		if ( ! is_set )
		{
			MyPrintfRed("Error: image %s not found in shift_blur_hack",input_search_images_filename);
			DEBUG_ABORT;
		}
		wxPrintf("For image %s using x,y shifts of %3.3f, %3.3f\n",input_search_images_filename,shift_hack_x[0],shift_hack_y[0]);


	}

	//
	remove_npix_from_edge = myroundint(particle_radius_angstroms / pixel_size);
//	wxPrintf("Removing %d pixels around the edge.\n", remove_npix_from_edge);

	Image input_image;
	Image padded_reference;
	Image input_reconstruction;
	Image template_reconstruction;
	Image current_projection;
	Image padded_projection;

	Image projection_filter;

	Image max_intensity_projection;

	Image best_psi;
	Image best_theta;
	Image best_phi;
	Image best_defocus;
	Image best_pixel_size;

	Image correlation_pixel_sum_image;
	Image correlation_pixel_sum_of_squares_image;

	Image temp_image;

	input_image.ReadSlice(&input_search_image_file, 1);

	// Resize input image to be factorizable by small numbers
	original_input_image_x = input_image.logical_x_dimension;
	original_input_image_y = input_image.logical_y_dimension;
	factorizable_x = input_image.logical_x_dimension;
	factorizable_y = input_image.logical_y_dimension;

	bool DO_FACTORIZATION = true;
	bool MUST_BE_POWER_OF_TWO = false; // Required for half-preicision xforms
	bool MUST_BE_FACTOR_OF_FOUR = true; // May be faster
	const int max_number_primes = 6;
	int primes[max_number_primes] = {2,3,5,7,9,13};
	float max_reduction_by_fraction_of_reference = 0.000001f; // FIXME the cpu version is crashing when the image is reduced, but not the GPU
	float max_increas_by_fraction_of_image = 0.1f;
	int max_padding = 0; // To restrict histogram calculation
	float histogram_padding_trim_rescale; // scale the counts to

	// for 5760 this will return
	// 5832 2     2     2     3     3     3     3     3     3 - this is ~ 10% faster than the previous solution BUT
	if (DO_FACTORIZATION)
	{
	for ( i = 0; i < max_number_primes; i++ )
	{

		factor_result_neg = ReturnClosestFactorizedLower(original_input_image_x, primes[i], true, MUST_BE_FACTOR_OF_FOUR);
		factor_result_pos = ReturnClosestFactorizedUpper(original_input_image_x, primes[i], true, MUST_BE_FACTOR_OF_FOUR);

//		wxPrintf("i, result, score = %i %i %g\n", i, factor_result, logf(float(abs(i) + 100)) * factor_result);
		if ( (float)(original_input_image_x - factor_result_neg) < (float)input_reconstruction_file.ReturnXSize() * max_reduction_by_fraction_of_reference)
		{
			factorizable_x = factor_result_neg;
			break;
		}
		if ((float)(-original_input_image_x + factor_result_pos) < (float)input_image.logical_x_dimension * max_increas_by_fraction_of_image)
		{
			factorizable_x = factor_result_pos;
			break;
		}

	}
	factor_score = FLT_MAX;
	for ( i = 0; i < max_number_primes; i++ )
	{

		factor_result_neg = ReturnClosestFactorizedLower(original_input_image_y, primes[i], true, MUST_BE_FACTOR_OF_FOUR);
		factor_result_pos = ReturnClosestFactorizedUpper(original_input_image_y, primes[i], true, MUST_BE_FACTOR_OF_FOUR);


//		wxPrintf("i, result, score = %i %i %g\n", i, factor_result, logf(float(abs(i) + 100)) * factor_result);
		if ( (float)(original_input_image_y - factor_result_neg) < (float)input_reconstruction_file.ReturnYSize() * max_reduction_by_fraction_of_reference)
		{
			factorizable_y = factor_result_neg;
			break;
		}
		if ((float)(-original_input_image_y + factor_result_pos) < (float)input_image.logical_y_dimension * max_increas_by_fraction_of_image)
		{
			factorizable_y = factor_result_pos;
			break;
		}

	}
	if (factorizable_x - original_input_image_x > max_padding) max_padding = factorizable_x - original_input_image_x;
	if (factorizable_y - original_input_image_y > max_padding) max_padding = factorizable_y - original_input_image_y;


	wxPrintf("old x, y; new x, y = %i %i %i %i\n", input_image.logical_x_dimension, input_image.logical_y_dimension, factorizable_x, factorizable_y);


	input_image.Resize(factorizable_x, factorizable_y, 1, input_image.ReturnAverageOfRealValuesOnEdges());
	if ( ! is_power_of_two(factorizable_x) && is_power_of_two(factorizable_y) )
	{
		// The speedup in the FFT for better factorization is also dependent on the dimension. The full transform (in cufft anyway) is faster if the best dimension is on X.
		// TODO figure out how to check the case where there is no factor of two, but one dimension is still faster. Probably getting around to writing an explicit planning tool would be useful.

		is_rotated_by_90 = true;
		input_image.Rotate2DInPlaceBy90Degrees(true);
	}

	}
	padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	max_intensity_projection.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	best_psi.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	best_theta.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	best_phi.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	best_defocus.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	best_pixel_size.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	correlation_pixel_sum_image.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	correlation_pixel_sum_of_squares_image.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
	double *correlation_pixel_sum = new double[input_image.real_memory_allocated];
	double *correlation_pixel_sum_of_squares = new double[input_image.real_memory_allocated];

	padded_reference.SetToConstant(0.0f);
	max_intensity_projection.SetToConstant(-FLT_MAX);
	best_psi.SetToConstant(0.0f);
	best_theta.SetToConstant(0.0f);
	best_phi.SetToConstant(0.0f);
	best_defocus.SetToConstant(0.0f);

	ZeroDoubleArray(correlation_pixel_sum, input_image.real_memory_allocated);
	ZeroDoubleArray(correlation_pixel_sum_of_squares, input_image.real_memory_allocated);

	input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices());
	if (padding != 1.0f)
	{
		input_reconstruction.Resize(input_reconstruction.logical_x_dimension * padding, input_reconstruction.logical_y_dimension * padding, input_reconstruction.logical_z_dimension * padding, input_reconstruction.ReturnAverageOfRealValuesOnEdges());
	}
//	input_reconstruction.ForwardFFT();
	//input_reconstruction.CosineMask(0.1, 0.01, true);
	//input_reconstruction.Whiten();
	//if (first_search_position == 0) input_reconstruction.QuickAndDirtyWriteSlices("/tmp/filter.mrc", 1, input_reconstruction.logical_z_dimension);
//	input_reconstruction.ZeroCentralPixel();
//	input_reconstruction.SwapRealSpaceQuadrants();

	sqrt_input_pixels =  sqrt((double)(input_image.logical_x_dimension * input_image.logical_y_dimension));
	// setup curve
	histogram_step = (histogram_max - histogram_min) / float(histogram_number_of_points);
	histogram_min_scaled = histogram_min / sqrt_input_pixels;
	histogram_step_scaled = histogram_step / sqrt_input_pixels;

	histogram_data = new long[histogram_number_of_points];

	for ( int counter = 0; counter < histogram_number_of_points; counter++ )
	{
		histogram_data[counter] = 0;
	}

	CTF input_ctf;
	input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));

	// assume cube

	current_projection.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);
	projection_filter.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);
	template_reconstruction.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_y_dimension, input_reconstruction.logical_z_dimension, true);
	if (padding != 1.0f) padded_projection.Allocate(input_reconstruction_file.ReturnXSize() * padding, input_reconstruction_file.ReturnXSize() * padding, false);


	// angular step
	float mask_radius_search;
	if (particle_radius_angstroms < 1.0f) { mask_radius_search = 200.0f; } // This was the original default value.
	else mask_radius_search = particle_radius_angstroms;

	if (angular_step <= 0)
	{
		angular_step = CalculateAngularStep(high_resolution_limit_search, mask_radius_search);
	}

	if (in_plane_angular_step <= 0)
	{
		psi_step = rad_2_deg(pixel_size / mask_radius_search);
		psi_step = 360.0 / int(360.0 / psi_step + 0.5);
	}
	else
	{
		psi_step = in_plane_angular_step;
	}

	//psi_start = psi_step / 2.0 * global_random_number_generator.GetUniformRandom();
	psi_start = 0.0f;
	psi_max = 360.0f;

	//psi_step = 5;

	//wxPrintf("psi_start = %f, psi_max = %f, psi_step = %f\n", psi_start, psi_max, psi_step);

	// search grid

	global_euler_search.InitGrid(my_symmetry, angular_step, 0.0f, 0.0f, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);
//	wxPrintf("%s",my_symmetry);
	if (my_symmetry.StartsWith("C1")) // TODO 2x check me - w/o this O symm at least is broken
	{
		if (global_euler_search.test_mirror == true) // otherwise the theta max is set to 90.0 and test_mirror is set to true.  However, I don't want to have to test the mirrors.
		{
			global_euler_search.theta_max = 180.0f;
		}
	}

	global_euler_search.CalculateGridSearchPositions(false);


	// for now, I am assuming the MTF has been applied already.
	// work out the filter to just whiten the image..

	whitening_filter.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
	number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

	wxDateTime my_time_out;
	wxDateTime my_time_in;

	// remove outliers
	input_image.ReplaceOutliersWithMean(13.0f); // This isn't valid for movie frames with poisson dist - can delete every non zero pixel
//	input_image.ReplaceOutliersWithMean(5.0f);
//	input_image.ReplacePoissonOutliersWithMode(5.0f);

	input_image.ForwardFFT();
	input_image.SwapRealSpaceQuadrants();

	input_image.ZeroCentralPixel();
	input_image.Compute1DPowerSpectrumCurve(&whitening_filter, &number_of_terms);
	whitening_filter.SquareRoot();
	whitening_filter.Reciprocal();
	whitening_filter.MultiplyByConstant(1.0f / whitening_filter.ReturnMaximumValue());

	//whitening_filter.WriteToFile("/tmp/filter.txt");
	input_image.ApplyCurveFilter(&whitening_filter);
	input_image.ZeroCentralPixel();
	input_image.DivideByConstant(sqrtf(input_image.ReturnSumOfSquares()));
	//input_image.QuickAndDirtyWriteSlice("/tmp/white.mrc", 1);
	//exit(-1);

	// count total searches (lazy)

	total_correlation_positions = 0;
	current_correlation_position = 0;

	// if running locally, search over all of them

	if (is_running_locally == true)
	{
		first_search_position = 0;
		last_search_position = global_euler_search.number_of_search_positions - 1;
	}

	// TODO unroll these loops and multiply the product.
	for (current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++)
	{
		//loop over each rotation angle

		for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
		{
			total_correlation_positions++;
		}
	}

	if (defocus_step <= 0.0)
	{
		defocus_search_range = 0.0f;
		defocus_step = 100.0f;
	}

	if (pixel_size_step <= 0.0f)
	{
		pixel_size_search_range = 0.0f;
		pixel_size_step = 0.02f;
	}

	total_correlation_positions *= (2 * myroundint(float(defocus_search_range)/float(defocus_step)) + 1);
	total_correlation_positions_per_thread = total_correlation_positions;

	number_of_rotations = 0;

	for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
	{
		number_of_rotations++;
	}

	ProgressBar *my_progress;

	//Loop over ever search position

	wxPrintf("\n\tFor image id %i\n",image_number_for_gui);
	wxPrintf("Searching %i positions on the Euler sphere (first-last: %i-%i)\n", last_search_position - first_search_position, first_search_position, last_search_position);
	wxPrintf("Searching %i rotations per position.\n", number_of_rotations);
	wxPrintf("There are %li correlation positions total.\n\n", total_correlation_positions);

	wxPrintf("Performing Search...\n\n");

//	wxPrintf("Searching %i - %i of %i total positions\n", first_search_position, last_search_position, global_euler_search.number_of_search_positions);
//	wxPrintf("psi_start = %f, psi_max = %f, psi_step = %f\n", psi_start, psi_max, psi_step);

	actual_number_of_ccs_calculated = 0.0;

	wxDateTime 	overall_start;
	wxDateTime 	overall_finish;
	overall_start = wxDateTime::Now();

	// These vars are only needed in the GPU code, but also need to be set out here to compile.
	bool first_gpu_loop = true;
	int nGPUs = 1;
	int nJobs = last_search_position-first_search_position+1;
	if (use_gpu && max_threads > nJobs)
	{
		wxPrintf("\n\tWarning, you request more threads (%d) than there are search positions (%d)\n", max_threads, nJobs);
		max_threads = nJobs;
	}

	int minPos = first_search_position;
	int maxPos = last_search_position;
	int incPos = (nJobs) / (max_threads);

//	wxPrintf("First last and inc %d, %d, %d\n", minPos, maxPos, incPos);
#ifdef ENABLEGPU
	TemplateMatchingCore *GPU;
	DeviceManager gpuDev;
#endif

	if (use_gpu)
	{
		total_correlation_positions_per_thread = total_correlation_positions / max_threads;

#ifdef ENABLEGPU
//	checkCudaErrors(cudaGetDeviceCount(&nGPUs));
	GPU = new TemplateMatchingCore[max_threads];
	gpuDev.Init(nGPUs);

//	wxPrintf("Host: %s is running\nnThreads: %d\nnGPUs: %d\n:nSearchPos %d \n",hostNameBuffer,nThreads, nGPUs, maxPos);

//	TemplateMatchingCore GPU(number_of_jobs_per_image_in_gui);
#endif
	}

	if (is_running_locally == true)
	{
		my_progress = new ProgressBar(total_correlation_positions_per_thread);
	}


//	wxPrintf("Starting job\n");
	for (size_i = - myroundint(float(pixel_size_search_range)/float(pixel_size_step)); size_i <= myroundint(float(pixel_size_search_range)/float(pixel_size_step)); size_i++)
	{


//		template_reconstruction.CopyFrom(&input_reconstruction);
		input_reconstruction.ChangePixelSize(&template_reconstruction, (pixel_size + float(size_i) * pixel_size_step) / pixel_size, 0.001f, true);
	//	template_reconstruction.ForwardFFT();
		template_reconstruction.ZeroCentralPixel();
		template_reconstruction.SwapRealSpaceQuadrants();

//		wxPrintf("First search last search position %d/ %d\n",first_search_position, last_search_position);

		if (use_gpu)
		{
#ifdef ENABLEGPU

	#pragma omp parallel num_threads(max_threads)
	{
		int tIDX = ReturnThreadNumberOfCurrentThread();
		gpuDev.SetGpu(tIDX);

		if (first_gpu_loop)
		{

				int t_first_search_position = first_search_position + (tIDX*incPos);
				int t_last_search_position = first_search_position + (incPos-1) + (tIDX*incPos);

				if (tIDX == (max_threads - 1)) t_last_search_position = maxPos;

				GPU[tIDX].Init(this, template_reconstruction, input_image, current_projection,
								pixel_size_search_range, pixel_size_step, pixel_size,
								defocus_search_range, defocus_step, defocus1, defocus2,
								psi_max, psi_start, psi_step,
								angles, global_euler_search,
								histogram_min_scaled, histogram_step_scaled,histogram_number_of_points,
								max_padding, t_first_search_position, t_last_search_position,
								my_progress, total_correlation_positions_per_thread, is_running_locally);

				wxPrintf("%d\n",tIDX);
				wxPrintf("%d\n", t_first_search_position);
				wxPrintf("%d\n", t_last_search_position);
				wxPrintf("Staring TemplateMatchingCore object %d to work on position range %d-%d\n", tIDX, t_first_search_position, t_last_search_position);

			first_gpu_loop = false;

		}
		else
		{
			GPU[tIDX].template_reconstruction.CopyFrom(&template_reconstruction);
		}
	} // end of omp block
#endif
	}
		for (defocus_i = - myroundint(float(defocus_search_range)/float(defocus_step)); defocus_i <= myroundint(float(defocus_search_range)/float(defocus_step)); defocus_i++)
		{


			// make the projection filter, which will be CTF * whitening filter
			input_ctf.SetDefocus((defocus1 + float(defocus_i) * defocus_step) / pixel_size, (defocus2 + float(defocus_i) * defocus_step) / pixel_size, deg_2_rad(defocus_angle));
//			input_ctf.SetDefocus((defocus1 + 200) / pixel_size, (defocus2 + 200) / pixel_size, deg_2_rad(defocus_angle));
			bool use_ctf_envelope = true;
			if (use_ctf_envelope)
			{
				wxPrintf("Using the CTF envelope!\n");
				input_ctf.SetEnvelope(voltage_kV, pixel_size, 8.0f / (pixel_size*pixel_size));
				projection_filter.CalculateCTFImage(input_ctf, false, true);

			}
			else projection_filter.CalculateCTFImage(input_ctf);
			projection_filter.ApplyCurveFilter(&whitening_filter);


			if (do_shift_blur_hack)
			{
				wxPrintf("Doing the shift blur hack loop\n");
				Image buffer_1;
				buffer_1.CopyFrom(&projection_filter);
				buffer_1.SetToConstant(0.0f);

				wxPrintf("Using a voltage and pixel size of %f %f\n",voltage_kV,pixel_size);
				ElectronDose my_electron_dose(voltage_kV, pixel_size);
				float *dose_filter = new float[projection_filter.real_memory_allocated/2];

//				angles.Init(global_euler_search.list_of_search_parameters[0][0], global_euler_search.list_of_search_parameters[0][1], 0, 0.0, 0.0);
////					angles.Init(130.0, 30.0, 199.5, 0.0, 0.0);
//				template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
//				current_projection.SwapRealSpaceQuadrants();
//
//				Image buffer_2 ;
//				Image buffer_3;
//				std::string fout;

				int i;
				int j;

				long pixel_counter = 0;

				float x_coordinate_2d;
				float y_coordinate_2d;


				float sinc_weight;

				for (int iFrame = 0; iFrame < blur_hack_n_frames; iFrame++)
				{

					ZeroFloatArray(dose_filter, projection_filter.real_memory_allocated/2);
					pixel_counter = 0;

					// Normally the pre-exposure is added to each frame. Here it is taken to be the total exposure.
					wxPrintf("for frame %d filtering to exposure %f\n",iFrame,shift_hack_d[iFrame]);
					my_electron_dose.CalculateDoseFilterAs1DArray(&projection_filter, dose_filter, 0.0f, shift_hack_d[iFrame]);

					for (j = 0; j <= projection_filter.physical_upper_bound_complex_y; j++)
					{

						y_coordinate_2d = projection_filter.ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * projection_filter.fourier_voxel_size_y / pixel_size;
						y_coordinate_2d *= shift_hack_y[iFrame];
						for (i = 0; i <= projection_filter.physical_upper_bound_complex_x; i++)
						{

							x_coordinate_2d = i * projection_filter.fourier_voxel_size_x / pixel_size;
							x_coordinate_2d *= shift_hack_x[iFrame];

							sinc_weight = sinc(1.0f*PIf*(x_coordinate_2d+y_coordinate_2d));
	//						if (sinc_weight < 0.f) sinc_weight *= sinc_weight; // make positive and shrink.
							if (do_exposure_filter_hack)
							{


								buffer_1.complex_values[pixel_counter] += (projection_filter.complex_values[pixel_counter] * sinc_weight * sqrtf(dose_filter[pixel_counter]));

							}
							else
							{
								buffer_1.complex_values[pixel_counter] += (projection_filter.complex_values[pixel_counter] * sinc_weight);

							}


							pixel_counter++;
						}
					}

//					buffer_2.CopyFrom(&current_projection);
//					buffer_3.CopyFrom(&buffer_1);
////					buffer_3.MultiplyByConstant(1./(1.+iFrame));
//					buffer_2.MultiplyPixelWise(buffer_3);
//					buffer_2.BackwardFFT();
//				    fout = "filtered_sinc_exp_" + std::to_string(iFrame) + ".mrc";
//					buffer_2.QuickAndDirtyWriteSlice(fout,1,false,1.5);

				}
//exit(-1);
				projection_filter.CopyFrom(&buffer_1);
				projection_filter.MultiplyByConstant(1.f/(float)blur_hack_n_frames);
//blur_hack_n_frames;

				delete [] dose_filter;

			}
//			projection_filter.QuickAndDirtyWriteSlice("normal_filter_with_envelope_and_sinc.mrc", 1, false, 1.5);
//			exit(0);
//			projection_filter.QuickAndDirtyWriteSlices("/tmp/projection_filter.mrc",1,projection_filter.logical_z_dimension,true,1.5);
			if (use_gpu)
			{
#ifdef ENABLEGPU
//			wxPrintf("\n\n\t\tsizeI defI %d %d\n\n\n", size_i, defocus_i);


			#pragma omp parallel num_threads(max_threads)
			{
				int tIDX = ReturnThreadNumberOfCurrentThread();
				gpuDev.SetGpu(tIDX);


				GPU[tIDX].RunInnerLoop(projection_filter, size_i, defocus_i, tIDX, current_correlation_position);



				#pragma omp critical
				{


					Image mip_buffer; mip_buffer.CopyFrom(&max_intensity_projection);
					Image psi_buffer; psi_buffer.CopyFrom(&max_intensity_projection);
					Image phi_buffer; phi_buffer.CopyFrom(&max_intensity_projection);
					Image theta_buffer; theta_buffer.CopyFrom(&max_intensity_projection);

					GPU[tIDX].d_max_intensity_projection.CopyDeviceToHost(mip_buffer, true, false);
					GPU[tIDX].d_best_psi.CopyDeviceToHost(psi_buffer, true, false);
					GPU[tIDX].d_best_phi.CopyDeviceToHost(phi_buffer, true, false);
					GPU[tIDX].d_best_theta.CopyDeviceToHost(theta_buffer, true, false);

//					mip_buffer.QuickAndDirtyWriteSlice("/tmp/tmpMipBuffer.mrc",1,1);
					// TODO should prob aggregate these across all workers
				// TODO add a copySum method that allocates a pinned buffer, copies there then sumes into the wanted image.
					Image sum;
					Image sumSq;

					sum.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
					sumSq.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);


					sum.SetToConstant(0.0f);
					sumSq.SetToConstant(0.0f);


					GPU[tIDX].d_sum3.CopyDeviceToHost(sum,true,false);
					GPU[tIDX].d_sumSq3.CopyDeviceToHost(sumSq,true,false);


					GPU[tIDX].d_max_intensity_projection.Wait();

					// TODO swap max_padding for explicit padding in x/y and limit calcs to that region.
					pixel_counter = 0;
					for (current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++)
					{
						for (current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++)
						{
							// first mip

							if (mip_buffer.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter])
							{
								max_intensity_projection.real_values[pixel_counter] = mip_buffer.real_values[pixel_counter];
								best_psi.real_values[pixel_counter] = psi_buffer.real_values[pixel_counter];
								best_theta.real_values[pixel_counter] = theta_buffer.real_values[pixel_counter];
								best_phi.real_values[pixel_counter] = phi_buffer.real_values[pixel_counter];
								best_defocus.real_values[pixel_counter] = float(defocus_i) * defocus_step;
								best_pixel_size.real_values[pixel_counter] = float(size_i) * pixel_size_step;

							}

							correlation_pixel_sum[pixel_counter] += (double)sum.real_values[pixel_counter];
							correlation_pixel_sum_of_squares[pixel_counter] += (double)sumSq.real_values[pixel_counter];

							pixel_counter++;
						}

						pixel_counter += max_intensity_projection.padding_jump_value;
					}


					GPU[tIDX].histogram.CopyToHostAndAdd(histogram_data);

//					current_correlation_position += GPU[tIDX].total_number_of_cccs_calculated;
					actual_number_of_ccs_calculated += GPU[tIDX].total_number_of_cccs_calculated;

				} // end of omp critical block
			} // end of parallel block


			continue;


#endif
			}

			for (current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++)
			{
				//loop over each rotation angle

				//current_rotation = 0;
				for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
				{

					angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], current_psi, 0.0, 0.0);
//					angles.Init(130.0, 30.0, 199.5, 0.0, 0.0);

					if (padding != 1.0f)
					{
						template_reconstruction.ExtractSlice(padded_projection, angles, 1.0f, false);
						padded_projection.SwapRealSpaceQuadrants();
						padded_projection.BackwardFFT();
						padded_projection.ClipInto(&current_projection);
						current_projection.ForwardFFT();
					}
					else
					{
						template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
						current_projection.SwapRealSpaceQuadrants();
					}
//					current_projection.QuickAndDirtyWriteSlice("proj.mrc", 1);
					//if (first_search_position == 0) current_projection.QuickAndDirtyWriteSlice("/tmp/small_proj_nofilter.mrc", 1);

					current_projection.MultiplyPixelWise(projection_filter);

					//if (first_search_position == 0) projection_filter.QuickAndDirtyWriteSlice("/tmp/projection_filter.mrc", 1);
					//if (first_search_position == 0) current_projection.QuickAndDirtyWriteSlice("/tmp/small_proj_afterfilter.mrc", 1);

					//current_projection.ZeroCentralPixel();
					//current_projection.DivideByConstant(sqrt(current_projection.ReturnSumOfSquares()));
					current_projection.BackwardFFT();
					//current_projection.ReplaceOutliersWithMean(6.0f);

					// find the pixel with the largest absolute density, and shift it to the centre

				/*	pixel_counter = 0;
					int best_x;
					int best_y;
					float max_value = -FLT_MAX;

					for ( int y = 0; y < current_projection.logical_y_dimension; y ++ )
					{
						for ( int x = 0; x < current_projection.logical_x_dimension; x ++ )
						{
							if (fabsf(current_projection.real_values[pixel_counter]) > max_value)
							{
								max_value = fabsf(current_projection.real_values[pixel_counter]);
								best_x = x - current_projection.physical_address_of_box_center_x;
								best_y = y - current_projection.physical_address_of_box_center_y;;
							}
							pixel_counter++;
						}
						pixel_counter += current_projection.padding_jump_value;
					}

					current_projection.RealSpaceIntegerShift(best_x, best_y, 0);
	*/
					///


					current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges());


//					variance = current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels;
//					current_projection.DivideByConstant(sqrtf(variance));
//					variance = current_projection.ReturnSumOfSquares();
					variance = current_projection.ReturnSumOfSquares() * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels \
							- powf(current_projection.ReturnAverageOfRealValues() * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
					current_projection.DivideByConstant(sqrtf(variance));
					current_projection.ClipIntoLargerRealSpace2D(&padded_reference);

					padded_reference.ForwardFFT();
					// Zeroing the central pixel is probably not doing anything useful...
					padded_reference.ZeroCentralPixel();
//					padded_reference.DivideByConstant(sqrtf(variance));

					//if (first_search_position == 0)  padded_reference.QuickAndDirtyWriteSlice("/tmp/proj.mrc", 1);

#ifdef MKL
					// Use the MKL
					vmcMulByConj(padded_reference.real_memory_allocated/2,reinterpret_cast <MKL_Complex8 *> (input_image.complex_values),reinterpret_cast <MKL_Complex8 *> (padded_reference.complex_values),reinterpret_cast <MKL_Complex8 *> (padded_reference.complex_values),VML_EP|VML_FTZDAZ_ON|VML_ERRMODE_IGNORE);
#else
					for (pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter ++)
					{
						padded_reference.complex_values[pixel_counter] = conj(padded_reference.complex_values[pixel_counter]) * input_image.complex_values[pixel_counter];
					}
#endif

					padded_reference.BackwardFFT();
//					padded_reference.QuickAndDirtyWriteSlice("cc.mrc", 1);
//					exit(0);

//					for (pixel_counter = 0; pixel_counter <  padded_reference.real_memory_allocated; pixel_counter++)
//					{
//						temp_float = padded_reference.real_values[pixel_counter] / variance;
//						padded_reference.real_values[pixel_counter] = temp_float * padded_reference.real_values[pixel_counter] - powf(temp_float, 2) * variance;
////						if (pixel_counter == 1000) wxPrintf("l, value = %g %g\n", temp_float, padded_reference.real_values[pixel_counter]);
////						padded_reference.real_values[pixel_counter] *= powf(float(padded_reference.number_of_real_space_pixels), 2);
////						padded_reference.real_values[pixel_counter] = temp_float;
//					}

					// update mip, and histogram..
					pixel_counter = 0;

					for (current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++)
					{
						for (current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++)
						{
							// first mip

							if (padded_reference.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter])
							{
								max_intensity_projection.real_values[pixel_counter] = padded_reference.real_values[pixel_counter];
								best_psi.real_values[pixel_counter] = current_psi;
								best_theta.real_values[pixel_counter] = global_euler_search.list_of_search_parameters[current_search_position][1];
								best_phi.real_values[pixel_counter] = global_euler_search.list_of_search_parameters[current_search_position][0];
								best_defocus.real_values[pixel_counter] = float(defocus_i) * defocus_step;
								best_pixel_size.real_values[pixel_counter] = float(size_i) * pixel_size_step;
//								if (size_i != 0) wxPrintf("size_i = %i\n", size_i);
//								correlation_pixel_sum[pixel_counter] = variance;
							}

							// histogram

							current_bin = int(double((padded_reference.real_values[pixel_counter]) - histogram_min_scaled) / histogram_step_scaled);
							//current_bin = int(double((padded_reference.real_values[pixel_counter]) - histogram_min) / histogram_step);

							if (current_bin >= 0 && current_bin <= histogram_number_of_points)
							{
								histogram_data[current_bin] += 1;
							}

							pixel_counter++;
						}

						pixel_counter+=padded_reference.padding_jump_value;
					}


//					correlation_pixel_sum.AddImage(&padded_reference);
					for (pixel_counter = 0; pixel_counter <  padded_reference.real_memory_allocated; pixel_counter++)
					{
						correlation_pixel_sum[pixel_counter] += padded_reference.real_values[pixel_counter];
					}
					padded_reference.SquareRealValues();
//					correlation_pixel_sum_of_squares.AddImage(&padded_reference);
					for (pixel_counter = 0; pixel_counter <  padded_reference.real_memory_allocated; pixel_counter++)
					{
						correlation_pixel_sum_of_squares[pixel_counter] += padded_reference.real_values[pixel_counter];
					}

					//max_intensity_projection.QuickAndDirtyWriteSlice("/tmp/mip.mrc", 1);

					current_projection.is_in_real_space = false;
					padded_reference.is_in_real_space = true;

					current_correlation_position++;
					if (is_running_locally == true) my_progress->Update(current_correlation_position);

					if (is_running_locally == false)
					{
						actual_number_of_ccs_calculated++;
						temp_float = current_correlation_position;
						JobResult *temp_result = new JobResult;
						temp_result->SetResult(1, &temp_float);
						AddJobToResultQueue(temp_result);
					}
				}
			}
		}
	}

	if (do_shift_blur_hack)
	{
		delete [] shift_hack_x;
		delete [] shift_hack_y;
		delete [] shift_hack_d;
	}

	wxPrintf("\n\n\tTimings: Overall: %s\n",(wxDateTime::Now()-overall_start).Format());



	for (pixel_counter = 0; pixel_counter <  input_image.real_memory_allocated; pixel_counter++)
	{
		correlation_pixel_sum_image.real_values[pixel_counter] = (float)correlation_pixel_sum[pixel_counter];
		correlation_pixel_sum_of_squares_image.real_values[pixel_counter] = (float)correlation_pixel_sum_of_squares[pixel_counter];
	}

	if (is_rotated_by_90)
	{
		// swap back all the images prior to re-sizing
		input_image.Rotate2DInPlaceBy90Degrees(false);
		max_intensity_projection.Rotate2DInPlaceBy90Degrees(false);

		best_psi.Rotate2DInPlaceBy90Degrees(false);
		best_theta.Rotate2DInPlaceBy90Degrees(false);
		best_phi.Rotate2DInPlaceBy90Degrees(false);
		best_defocus.Rotate2DInPlaceBy90Degrees(false);
		best_pixel_size.Rotate2DInPlaceBy90Degrees(false);

		correlation_pixel_sum_image.Rotate2DInPlaceBy90Degrees(false);
		correlation_pixel_sum_of_squares_image.Rotate2DInPlaceBy90Degrees(false);

		// This is ineffecient, but a quick way to ensure consistent results.
		delete [] correlation_pixel_sum;
		delete [] correlation_pixel_sum_of_squares;
		// Now we have the rotated values which may also be a different total amount of memory
		correlation_pixel_sum = new double[input_image.real_memory_allocated];
		correlation_pixel_sum_of_squares = new double[input_image.real_memory_allocated];
		ZeroDoubleArray(correlation_pixel_sum, input_image.real_memory_allocated);
		ZeroDoubleArray(correlation_pixel_sum_of_squares, input_image.real_memory_allocated);
		for (pixel_counter = 0; pixel_counter <  input_image.real_memory_allocated; pixel_counter++)
		{
			correlation_pixel_sum[pixel_counter] = (double)correlation_pixel_sum_image.real_values[pixel_counter];
			correlation_pixel_sum_of_squares[pixel_counter] = (double)correlation_pixel_sum_of_squares_image.real_values[pixel_counter];
		}

	}



	if (is_running_locally == true)
	{
		delete my_progress;

		// scale images..

		for (pixel_counter = 0; pixel_counter <  input_image.real_memory_allocated; pixel_counter++)
		{

//			correlation_pixel_sum.real_values[pixel_counter] /= float(total_correlation_positions);
//			correlation_pixel_sum_of_squares.real_values[pixel_counter] = correlation_pixel_sum_of_squares.real_values[pixel_counter] / float(total_correlation_positions) - powf(correlation_pixel_sum.real_values[pixel_counter], 2);
//			if (correlation_pixel_sum_of_squares.real_values[pixel_counter] > 0.0f)
//			{
//				correlation_pixel_sum_of_squares.real_values[pixel_counter] = sqrtf(correlation_pixel_sum_of_squares.real_values[pixel_counter]) * sqrtf(correlation_pixel_sum.logical_x_dimension * correlation_pixel_sum.logical_y_dimension);
//			}
//			else correlation_pixel_sum_of_squares.real_values[pixel_counter] = 0.0f;
			correlation_pixel_sum[pixel_counter] /= float(total_correlation_positions);
			correlation_pixel_sum_of_squares[pixel_counter] = correlation_pixel_sum_of_squares[pixel_counter] / float(total_correlation_positions) - powf(correlation_pixel_sum[pixel_counter], 2);
			if (correlation_pixel_sum_of_squares[pixel_counter] > 0.0f)
			{
				correlation_pixel_sum_of_squares[pixel_counter] = sqrtf(correlation_pixel_sum_of_squares[pixel_counter]) * (float)sqrt_input_pixels;
			}
			else correlation_pixel_sum_of_squares[pixel_counter] = 0.0f;
			correlation_pixel_sum[pixel_counter] *= (float)sqrt_input_pixels;

		}


		max_intensity_projection.MultiplyByConstant((float)sqrt_input_pixels);
//		correlation_pixel_sum.MultiplyByConstant(sqrtf(max_intensity_projection.logical_x_dimension * max_intensity_projection.logical_y_dimension));
//		correlation_pixel_sum_of_squares.MultiplyByConstant(max_intensity_projection.logical_x_dimension * max_intensity_projection.logical_y_dimension);

		// we need to quadrant swap the images, also shift them, with an extra pixel shift.  This is because I take the conjugate of the input image, not the reference..



//		max_intensity_projection.InvertPixelOrder();
//		max_intensity_projection.SwapRealSpaceQuadrants();


//		best_psi.InvertPixelOrder();
//		best_psi.SwapRealSpaceQuadrants();

//		best_theta.InvertPixelOrder();
//		best_theta.SwapRealSpaceQuadrants();

//		best_phi.InvertPixelOrder();
//		best_phi.SwapRealSpaceQuadrants();

//		best_defocus.InvertPixelOrder();
//		best_defocus.SwapRealSpaceQuadrants();

//		correlation_pixel_sum.InvertPixelOrder();
//		correlation_pixel_sum.SwapRealSpaceQuadrants();

//		correlation_pixel_sum_of_squares.InvertPixelOrder();
//		correlation_pixel_sum_of_squares.SwapRealSpaceQuadrants();



		// calculate the expected threshold (from peter's paper)
		const float CCG_NOISE_STDDEV = 1.0;
		double temp_threshold;
		double erf_input = 2.0 / (1.0 * (double)original_input_image_x * (double)original_input_image_y * (double)total_correlation_positions);
#ifdef MKL
		vdErfcInv(1, &erf_input, &temp_threshold);
#else
		cisTEM_erfcinv(erf_input);
#endif
		expected_threshold = sqrtf(2.0f)*(float)temp_threshold*CCG_NOISE_STDDEV;

//		expected_threshold = sqrtf(2.0f)*cisTEM_erfcinv((2.0f*(1))/((original_input_image_x * original_input_image_y * double(total_correlation_positions))));




		// write out images..

//		wxPrintf("\nPeak at %g, %g : %g\n", max_intensity_projection.FindPeakWithIntegerCoordinates().x, max_intensity_projection.FindPeakWithIntegerCoordinates().y, max_intensity_projection.FindPeakWithIntegerCoordinates().value);
//		wxPrintf("Sigma = %g, ratio = %g\n", sqrtf(max_intensity_projection.ReturnVarianceOfRealValues()), max_intensity_projection.FindPeakWithIntegerCoordinates().value / sqrtf(max_intensity_projection.ReturnVarianceOfRealValues()));

		temp_image.CopyFrom(&max_intensity_projection);
		temp_image.Resize(original_input_image_x, original_input_image_y, 1, temp_image.ReturnAverageOfRealValuesOnEdges());
		temp_image.QuickAndDirtyWriteSlice(mip_output_file.ToStdString(), 1, pixel_size);
//		max_intensity_projection.SubtractImage(&correlation_pixel_sum);
		for (pixel_counter = 0; pixel_counter <  input_image.real_memory_allocated; pixel_counter++)
		{
			max_intensity_projection.real_values[pixel_counter] -= correlation_pixel_sum[pixel_counter];
			if (correlation_pixel_sum_of_squares[pixel_counter] > 0.0f)
			{
				max_intensity_projection.real_values[pixel_counter] /= correlation_pixel_sum_of_squares[pixel_counter];
			}
			else max_intensity_projection.real_values[pixel_counter] = 0.0f;
			correlation_pixel_sum_image.real_values[pixel_counter] = correlation_pixel_sum[pixel_counter];
			correlation_pixel_sum_of_squares_image.real_values[pixel_counter] = correlation_pixel_sum_of_squares[pixel_counter];
		}
//		max_intensity_projection.DividePixelWise(correlation_pixel_sum_of_squares);
		max_intensity_projection.Resize(original_input_image_x, original_input_image_y, 1, max_intensity_projection.ReturnAverageOfRealValuesOnEdges());
		max_intensity_projection.QuickAndDirtyWriteSlice(scaled_mip_output_file.ToStdString(), 1, pixel_size);


		correlation_pixel_sum_image.Resize(original_input_image_x, original_input_image_y, 1, correlation_pixel_sum_image.ReturnAverageOfRealValuesOnEdges());
		correlation_pixel_sum_image.QuickAndDirtyWriteSlice(correlation_avg_output_file.ToStdString(), 1, pixel_size);
		correlation_pixel_sum_of_squares_image.Resize(original_input_image_x, original_input_image_y, 1, correlation_pixel_sum_of_squares_image.ReturnAverageOfRealValuesOnEdges());
		correlation_pixel_sum_of_squares_image.QuickAndDirtyWriteSlice(correlation_std_output_file.ToStdString(), 1, pixel_size);
		best_psi.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
		best_psi.QuickAndDirtyWriteSlice(best_psi_output_file.ToStdString(), 1, pixel_size);
		best_theta.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
		best_theta.QuickAndDirtyWriteSlice(best_theta_output_file.ToStdString(), 1, pixel_size);
		best_phi.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
		best_phi.QuickAndDirtyWriteSlice(best_phi_output_file.ToStdString(), 1, pixel_size);
		best_defocus.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
		best_defocus.QuickAndDirtyWriteSlice(best_defocus_output_file.ToStdString(), 1, pixel_size);
		best_pixel_size.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
		best_pixel_size.QuickAndDirtyWriteSlice(best_pixel_size_output_file.ToStdString(), 1, pixel_size);

		// write out histogram..

		temp_float = histogram_min + (histogram_step / 2.0f); // start position
		NumericTextFile histogram_file(output_histogram_file, OPEN_TO_WRITE, 4);

		double *expected_survival_histogram = new double[histogram_number_of_points];
		double *survival_histogram = new double[histogram_number_of_points];
		ZeroDoubleArray(survival_histogram, histogram_number_of_points);

		for (int line_counter = 0; line_counter <= histogram_number_of_points; line_counter++)
		{
				expected_survival_histogram[line_counter] = (erfc((temp_float + histogram_step * float(line_counter))/sqrtf(2.0f))/2.0f)*((float)(sqrt_input_pixels*sqrt_input_pixels) * float(total_correlation_positions));
		}

		survival_histogram[histogram_number_of_points - 1] = histogram_data[histogram_number_of_points - 1];

		for (int line_counter = histogram_number_of_points - 2; line_counter >= 0 ; line_counter--)
		{
			survival_histogram[line_counter] = survival_histogram[line_counter + 1] + histogram_data[line_counter];
		}

		histogram_file.WriteCommentLine("Expected threshold = %.2f\n", expected_threshold);
		histogram_file.WriteCommentLine("SNR, histogram, survival histogram, random survival histogram");

		for (int line_counter = 0; line_counter < histogram_number_of_points; line_counter++)
		{
			temp_double_array[0] = temp_float + histogram_step * float(line_counter);
			temp_double_array[1] = histogram_data[line_counter];
			temp_double_array[2] = survival_histogram[line_counter];
			temp_double_array[3] = expected_survival_histogram[line_counter];
			histogram_file.WriteLine(temp_double_array);
		}

		histogram_file.Close();

		// memory cleanup

		delete [] survival_histogram;
		delete [] expected_survival_histogram;
	}
	else
	{
		// send back the final images to master (who should merge them, and send to the gui)

		long result_array_counter;
		long number_of_result_floats = number_of_meta_data_values; // first float is x size, 2nd is y size of images, 3rd is number allocated, 4th  float is number of doubles in the histogram
		long pixel_counter;
		float *pointer_to_histogram_data;

		pointer_to_histogram_data = (float *) histogram_data;

//		max_intensity_projection.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		correlation_pixel_sum_image.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		correlation_pixel_sum_of_squares_image.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		best_psi.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		best_theta.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		best_phi.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		best_defocus.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);
//		best_pixel_size.Resize(original_input_image_x, original_input_image_y, 1, 0.0f);

		// If the padded image is large, we want to resize, then trim to valid area, otherwise we want to trim to valid area and then resize.
		// Default to the case where the padding increases the image size. A call to resize a same size image only cost the function call time.
		int	trim_x = original_input_image_x - remove_npix_from_edge;
		int	trim_y = original_input_image_y - remove_npix_from_edge;

		float central_average;
		float central_region = 0.35*(max_intensity_projection.logical_x_dimension + max_intensity_projection.logical_y_dimension - 2 * remove_npix_from_edge);

		if (original_input_image_x > max_intensity_projection.logical_x_dimension)
		{
			trim_x = max_intensity_projection.logical_x_dimension - remove_npix_from_edge;
		}
		if (original_input_image_y > max_intensity_projection.logical_y_dimension)
		{
			trim_y = max_intensity_projection.logical_y_dimension - remove_npix_from_edge;
		}


		// mip
		central_average = max_intensity_projection.ReturnAverageOfRealValues(central_region, false);
		max_intensity_projection.Resize(trim_x, trim_y, 1, central_average);

		max_intensity_projection.Resize(original_input_image_x, original_input_image_y, 1, central_average);


		//sum
		central_average = correlation_pixel_sum_image.ReturnAverageOfRealValues(central_region, false);
		correlation_pixel_sum_image.Resize(trim_x, trim_y, 1, central_average);
		correlation_pixel_sum_image.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// sq sum
		central_average = correlation_pixel_sum_of_squares_image.ReturnAverageOfRealValues(central_region, false);
		correlation_pixel_sum_of_squares_image.Resize(trim_x, trim_y, 1, central_average);
		correlation_pixel_sum_of_squares_image.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// psi
		central_average = best_psi.ReturnAverageOfRealValues(central_region, false);
		best_psi.Resize(trim_x, trim_y, 1, central_average);
		best_psi.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// theta
		central_average = best_theta.ReturnAverageOfRealValues(central_region, false);
		best_theta.Resize(trim_x, trim_y, 1, central_average);
		best_theta.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// phi
		central_average = best_phi.ReturnAverageOfRealValues(central_region, false);
		best_phi.Resize(trim_x, trim_y, 1, central_average);
		best_phi.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// pixel
		central_average = best_pixel_size.ReturnAverageOfRealValues(central_region, false);
		best_pixel_size.Resize(trim_x, trim_y, 1, central_average);
		best_pixel_size.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// defocus
		central_average = best_defocus.ReturnAverageOfRealValues(central_region, false);
		best_defocus.Resize(trim_x, trim_y, 1, central_average);
		best_defocus.Resize(original_input_image_x, original_input_image_y, 1, central_average);

		// Make sure there is enough space allocated for all results
		number_of_result_floats += max_intensity_projection.real_memory_allocated * number_of_output_images;
		number_of_result_floats += histogram_number_of_points * sizeof(long)/sizeof(float); // histogram are longs

		float *result = new float[number_of_result_floats];
		// Not zero floating this array since all additions are assignments. This can help to expose any indexing errors.

		result[0] = max_intensity_projection.logical_x_dimension;
		result[1] = max_intensity_projection.logical_y_dimension;
		result[2] = max_intensity_projection.real_memory_allocated;
		result[3] = histogram_number_of_points;
		result[4] = actual_number_of_ccs_calculated;
		result[5] = (float)sqrt_input_pixels;
//		result[5] = original_input_image_x;
//		result[6] = original_input_image_y;

		result_array_counter = number_of_meta_data_values;

		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = max_intensity_projection.real_values[pixel_counter];
			result_array_counter++;
		}

		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = best_psi.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = best_theta.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = best_phi.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = best_defocus.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = best_pixel_size.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = correlation_pixel_sum_image.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < max_intensity_projection.real_memory_allocated; pixel_counter++)
		{
			result[result_array_counter] = correlation_pixel_sum_of_squares_image.real_values[pixel_counter];
			result_array_counter++;
		}


		for (pixel_counter = 0; pixel_counter < histogram_number_of_points * 2; pixel_counter++)
		{
			result[result_array_counter] = 	pointer_to_histogram_data[pixel_counter];
			result_array_counter++;
		}

		SendProgramDefinedResultToMaster(result, number_of_result_floats, image_number_for_gui, number_of_jobs_per_image_in_gui);
	}

	delete [] histogram_data;

	if (is_running_locally == true)
	{
		wxPrintf("\nMatch Template: Normal termination\n");
		wxDateTime finish_time = wxDateTime::Now();
		wxPrintf("Total Run Time : %s\n\n", finish_time.Subtract(start_time).Format("%Hh:%Mm:%Ss"));
	}

	return true;
}

void MatchTemplateApp::MasterHandleProgramDefinedResult(float *result_array, long array_size, int result_number, int number_of_expected_results)
{
	// do we have this image number already?

	bool need_a_new_result = true;
	int array_location = -1;
	long pixel_counter;

	wxPrintf("Master Handling result for image %i..", result_number);

	for (int result_counter = 0; result_counter < aggregated_results.GetCount(); result_counter++)
	{
		if (aggregated_results[result_counter].image_number == result_number)
		{
			aggregated_results[result_counter].AddResult(result_array, array_size, result_number, number_of_expected_results);
			need_a_new_result = false;
			array_location = result_counter;
			wxPrintf("Found array location for image %i, at %i\n", result_number, array_location);
			break;
		}
	}

	if (need_a_new_result == true) // we aren't collecting data for this result yet.. start
	{
		AggregatedTemplateResult result_to_add;
		aggregated_results.Add(result_to_add);
		aggregated_results[aggregated_results.GetCount() - 1].image_number = result_number;
		aggregated_results[aggregated_results.GetCount() - 1].AddResult(result_array, array_size, result_number, number_of_expected_results);
		array_location = aggregated_results.GetCount() - 1;
		wxPrintf("Adding new result to array for image %i, at %i\n", result_number, array_location);
	}

	// did this complete a result?

	if (aggregated_results[array_location].number_of_received_results == number_of_expected_results) // we should be done for this image
	{
		// TODO send the result back to the GUI, for now hack mode to save the files to the directory..

		wxString directory_for_writing_results = current_job_package.jobs[0].arguments[37].ReturnStringArgument();

//		wxPrintf("temp x, y, n, resize x, y = %i %i %i %i %i \n", int(aggregated_results[array_location].collated_data_array[0]), \
//			int(aggregated_results[array_location].collated_data_array[1]), int(result_array[2]), int(result_array[5]), int(result_array[6]));

		Image temp_image;

		Image scaled_mip;
		Image psi_image;
		Image phi_image;
		Image theta_image;
		Image defocus_image;
		Image pixel_size_image;

		Image result_image;
		Image input_reconstruction;
		Image current_projection;

		int number_of_peaks_found = 0;
		float sq_dist_x;
		float sq_dist_y;
		float current_phi;
		float current_psi;
		float current_theta;
		int i;
		int j;
		long address;

		ArrayOfTemplateMatchFoundPeakInfos all_peak_infos;
		TemplateMatchFoundPeakInfo temp_peak_info;

		Peak current_peak;
		AnglesAndShifts angles;

		double sqrt_input_pixels = aggregated_results[array_location].collated_data_array[5];
		bool use_gpu = current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[40].ReturnBoolArgument();

		ImageFile input_reconstruction_file;
		input_reconstruction_file.OpenFile(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[1].ReturnStringArgument(), false);

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);

		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_mip_data[pixel_counter] * sqrt_input_pixels;
		}


		wxPrintf("Writing result %i\n", aggregated_results[array_location].image_number - 1);
		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[21].ReturnStringArgument(), 1);
		temp_image.Deallocate();

		// psi

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_psi_data[pixel_counter];
		}



		//temp_image.QuickAndDirtyWriteSlice(wxString::Format("%s/psi.mrc", directory_for_writing_results).ToStdString(), aggregated_results[array_location].image_number);
		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[22].ReturnStringArgument(), 1);
		psi_image.CopyFrom(&temp_image);
		temp_image.Deallocate();

		//theta

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_theta_data[pixel_counter];
		}


		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[23].ReturnStringArgument(), 1);
		theta_image.CopyFrom(&temp_image);
		temp_image.Deallocate();


		// phi

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_phi_data[pixel_counter];
		}


		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[24].ReturnStringArgument(), 1);
		phi_image.CopyFrom(&temp_image);
		temp_image.Deallocate();


		// defocus

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_defocus_data[pixel_counter];
		}


		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[25].ReturnStringArgument(), 1);
		defocus_image.CopyFrom(&temp_image);
		temp_image.Deallocate();


		// pixel size

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_pixel_size_data[pixel_counter];
		}


		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[26].ReturnStringArgument(), 1);
		pixel_size_image.CopyFrom(&temp_image);
		temp_image.Deallocate();

		// do the scaling...

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			aggregated_results[array_location].collated_pixel_sums[pixel_counter] /= aggregated_results[array_location].total_number_of_ccs;
			aggregated_results[array_location].collated_pixel_square_sums[pixel_counter] = sqrtf(aggregated_results[array_location].collated_pixel_square_sums[pixel_counter] /
																					       aggregated_results[array_location].total_number_of_ccs - powf(aggregated_results[array_location].collated_pixel_sums[pixel_counter], 2));
			if (aggregated_results[array_location].collated_pixel_square_sums[pixel_counter] > 0.0f)
			{

				// Save the variance, not the stdDev
//				aggregated_results[array_location].collated_pixel_square_sums[pixel_counter] = sqrtf(aggregated_results[array_location].collated_pixel_square_sums[pixel_counter]);
				aggregated_results[array_location].collated_mip_data[pixel_counter] = (aggregated_results[array_location].collated_mip_data[pixel_counter] - aggregated_results[array_location].collated_pixel_sums[pixel_counter]) /
																					   aggregated_results[array_location].collated_pixel_square_sums[pixel_counter];
			}
			else
			{
				aggregated_results[array_location].collated_pixel_square_sums[pixel_counter] = 0.0f;
				aggregated_results[array_location].collated_mip_data[pixel_counter] = 0.0f;
			}
		}

		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_mip_data[pixel_counter];
		}


		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[27].ReturnStringArgument(), 1);
		scaled_mip.CopyFrom(&temp_image);
		temp_image.Deallocate();


		// sums

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_pixel_sums[pixel_counter] * sqrt_input_pixels;
		}

		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[28].ReturnStringArgument(), 1);
		temp_image.Deallocate();


		// square sums

		temp_image.Allocate(int(aggregated_results[array_location].collated_data_array[0]), int(aggregated_results[array_location].collated_data_array[1]), true);
		for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
		{
			temp_image.real_values[pixel_counter] = aggregated_results[array_location].collated_pixel_square_sums[pixel_counter] * sqrt_input_pixels;
		}

		temp_image.QuickAndDirtyWriteSlice(current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[36].ReturnStringArgument(), 1);
		temp_image.Deallocate();


		// histogram

		float histogram_step = (histogram_max - histogram_min) / float(histogram_number_of_points);
		float temp_float = histogram_min + (histogram_step / 2.0f); // start position
		//NumericTextFile histogram_file(wxString::Format("%s/histogram_%i.txt", directory_for_writing_results, aggregated_results[array_location].image_number), OPEN_TO_WRITE, 4);
		NumericTextFile histogram_file(	current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[31].ReturnStringArgument(), OPEN_TO_WRITE, 4);

		double *expected_survival_histogram = new double[histogram_number_of_points];
		double *survival_histogram = new double[histogram_number_of_points];

		float expected_threshold;

		double temp_double_array[5];

		ZeroDoubleArray(survival_histogram, histogram_number_of_points);
		survival_histogram[histogram_number_of_points - 1] = aggregated_results[array_location].collated_histogram_data[histogram_number_of_points - 1];

		for (int line_counter = histogram_number_of_points - 2; line_counter >= 0 ; line_counter--)
		{
			survival_histogram[line_counter] = survival_histogram[line_counter + 1] + aggregated_results[array_location].collated_histogram_data[line_counter];
		}

		for (int line_counter = 0; line_counter <= histogram_number_of_points; line_counter++)
		{
			expected_survival_histogram[line_counter] = (erfc((temp_float + histogram_step * float(line_counter))/sqrtf(2.0f))/2.0f)*(aggregated_results[array_location].collated_data_array[0] * aggregated_results[array_location].collated_data_array[1] * aggregated_results[array_location].total_number_of_ccs);
		}

		// calculate the expected threshold (from peter's paper)
		const float CCG_NOISE_STDDEV = 1.0;
		double temp_threshold = 0.0;
		double erf_input = 2.0 / (1.0 * ((double)aggregated_results[array_location].collated_data_array[0] * (double)aggregated_results[array_location].collated_data_array[1] * (double)aggregated_results[array_location].total_number_of_ccs));
//		wxPrintf("ox oy total %3.3e %3.3e %3.3e\n", (double)result_array[5] , (double)result_array[6] , (double)aggregated_results[array_location].total_number_of_ccs, erf_input);

#ifdef MKL
		vdErfcInv(1, &erf_input, &temp_threshold);
#else
		temp_threshold = cisTEM_erfcinv(erf_input);
#endif
		expected_threshold = sqrtf(2.0f)*(float)temp_threshold*CCG_NOISE_STDDEV;

//		expected_threshold = sqrtf(2.0f)*cisTEM_erfcinv((2.0f*(1))/(((original_input_image_x * original_input_image_y * aggregated_results[array_location].total_number_of_ccs))));

		histogram_file.WriteCommentLine("Expected threshold = %.2f\n", expected_threshold);
		histogram_file.WriteCommentLine("histogram, expected histogram, survival histogram, expected survival histogram");

		if (use_gpu)
		{
		// In the GPU code, I am not histogramming the padding regions which are not valid. Adjust the counts here. Maybe not the best approach. FIXME also the cpu counts.
#ifdef ENABLEGPU
		double sum_expected = 0.0;
		double sum_counted = 0.0;

		for (int line_counter = 0; line_counter < histogram_number_of_points; line_counter++)
		{
			sum_counted += survival_histogram[line_counter];
			sum_expected += expected_survival_histogram[line_counter];
		}
		for (int line_counter = 0; line_counter < histogram_number_of_points; line_counter++)
		{
			survival_histogram[line_counter] *= (float)(sum_expected / sum_counted);
		}
#endif
		}


		for (int line_counter = 0; line_counter < histogram_number_of_points; line_counter++)
		{
			temp_double_array[0] = temp_float + histogram_step * float(line_counter);
			temp_double_array[1] = aggregated_results[array_location].collated_histogram_data[line_counter];
			temp_double_array[2] = survival_histogram[line_counter];
			temp_double_array[3] = expected_survival_histogram[line_counter];
			histogram_file.WriteLine(temp_double_array);
		}

		histogram_file.Close();

		// Calculate the result image, and keep the peak info to send back...

		int   min_peak_radius = current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[39].ReturnFloatArgument();
		float min_peak_radius_squared = powf(float(min_peak_radius), 2);


		result_image.Allocate(scaled_mip.logical_x_dimension, scaled_mip.logical_y_dimension, 1);
		result_image.SetToConstant(0.0f);

		input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices());
		float max_density = input_reconstruction.ReturnAverageOfMaxN();
		input_reconstruction.DivideByConstant(max_density);
		input_reconstruction.ForwardFFT();
		input_reconstruction.MultiplyByConstant(sqrtf(input_reconstruction.logical_x_dimension * input_reconstruction.logical_y_dimension * sqrtf(input_reconstruction.logical_z_dimension)));
		input_reconstruction.ZeroCentralPixel();
		input_reconstruction.SwapRealSpaceQuadrants();

		// assume cube

		current_projection.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_x_dimension, false);

		// loop until the found peak is below the threshold

		long nTrys = 0;
		while (1==1)
		{
			// look for a peak..
			nTrys++;
//			wxPrintf("Trying the %ld'th peak\n",nTrys);
			// FIXME min-distance from edges would be better to set dynamically.
			current_peak = scaled_mip.FindPeakWithIntegerCoordinates(0.0, FLT_MAX, input_reconstruction.logical_x_dimension / 4 + 1);
			if (current_peak.value < expected_threshold) break;

			// ok we have peak..

			number_of_peaks_found++;

			// get angles and mask out the local area so it won't be picked again..

			address = 0;

			current_peak.x = current_peak.x + scaled_mip.physical_address_of_box_center_x;
			current_peak.y = current_peak.y + scaled_mip.physical_address_of_box_center_y;

			// arguments[2] = pixel_size
			temp_peak_info.x_pos = current_peak.x * current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[2].ReturnFloatArgument(); // RETURNING IN ANGSTROMS
			temp_peak_info.y_pos = current_peak.y * current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[2].ReturnFloatArgument(); // RETURNING IN ANGSTROMS


//			wxPrintf("Peak = %f, %f, %f : %f\n", current_peak.x, current_peak.y, current_peak.value);

			for ( j = std::max(myroundint(current_peak.y) - min_peak_radius, 0); j < std::min(myroundint(current_peak.y) + min_peak_radius, scaled_mip.logical_y_dimension) ; j ++ )
			{
				sq_dist_y = float(j)-current_peak.y;
				sq_dist_y *= sq_dist_y;

				for ( i = std::max(myroundint(current_peak.x) - min_peak_radius, 0); i < std::min(myroundint(current_peak.x) + min_peak_radius, scaled_mip.logical_x_dimension) ; i ++ )
				{
					sq_dist_x = float(i)-current_peak.x;
					sq_dist_x *= sq_dist_x;
					address = phi_image.ReturnReal1DAddressFromPhysicalCoord(i,j,0);

					// The square centered at the pixel
					if (sq_dist_x == 0 && sq_dist_y == 0)
					{
						current_phi = phi_image.real_values[address];
						current_theta = theta_image.real_values[address];
						current_psi = psi_image.real_values[address];

						temp_peak_info.phi = phi_image.real_values[address];
						temp_peak_info.theta = theta_image.real_values[address];
						temp_peak_info.psi = psi_image.real_values[address];

						temp_peak_info.defocus = defocus_image.real_values[address];  // RETURNING MINUS
						temp_peak_info.pixel_size = pixel_size_image.real_values[address];
						temp_peak_info.peak_height = scaled_mip.real_values[address];

					}

					if ( sq_dist_x + sq_dist_y <= min_peak_radius_squared )
					{
						scaled_mip.real_values[address] = -FLT_MAX;
					}


//					address++;
				}
//				address += scaled_mip.padding_jump_value;
			}


	//		wxPrintf("Peak %4i at x, y, psi, theta, phi, defocus, pixel size = %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f : %10.6f\n", number_of_peaks_found, current_peak.x, current_peak.y, current_psi, current_theta, current_phi, current_defocus, current_pixel_size, current_peak.value);
	//		coordinates[0] = current_peak.x * pixel_size;
	//		coordinates[1] = current_peak.y * pixel_size;
	////		coordinates[2] = binned_pixel_size * (slab.physical_address_of_box_center_z - binned_reconstruction.physical_address_of_box_center_z) - current_defocus;
	//		coordinates[2] = binned_pixel_size * slab.physical_address_of_box_center_z - current_defocus;
	//		coordinate_file.WriteLine(coordinates);

			// ok get a projection

			//////////////////////////////////////////////
			// CURRENTLY HARD CODED TO ONLY DO 1000 MAX //
			//////////////////////////////////////////////

			if (number_of_peaks_found <= MAX_ALLOWED_NUMBER_OF_PEAKS)
			{

				angles.Init(current_phi, current_theta, current_psi, 0.0, 0.0);

				input_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
				current_projection.SwapRealSpaceQuadrants();

				current_projection.MultiplyByConstant(sqrtf(current_projection.logical_x_dimension * current_projection.logical_y_dimension));
				current_projection.BackwardFFT();
				current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges());

				// insert it into the output image

				result_image.InsertOtherImageAtSpecifiedPosition(&current_projection, current_peak.x - result_image.physical_address_of_box_center_x, current_peak.y - result_image.physical_address_of_box_center_y, 0, 0.0f);
				all_peak_infos.Add(temp_peak_info);

				//current_projection.QuickAndDirtyWriteSlice("/tmp/projs.mrc", all_peak_infos.GetCount());
			}
			else
			{
//				SendError("More than 1000 peaks above threshold were found. Limiting results to 1000 peaks.\n");
//				break;
				wxPrintf("Something seems to have gone wrong, more than 1000 _peaks_ were found\n");
				scaled_mip.QuickAndDirtyWriteSlice("/tmp/scaled_mip_1000.mrc", 1);
				exit(0);
			}

		}

		// save the output image

		result_image.QuickAndDirtyWriteSlice(	current_job_package.jobs[(aggregated_results[array_location].image_number - 1) * number_of_expected_results].arguments[38].ReturnStringArgument(), 1, true);

		// tell the gui that this result is available...

		ArrayOfTemplateMatchFoundPeakInfos blank_changes;
		SendTemplateMatchingResultToSocket(controller_socket, aggregated_results[array_location].image_number, expected_threshold, all_peak_infos, blank_changes);

				// this should be done now.. so delete it

		aggregated_results.RemoveAt(array_location);
		delete [] expected_survival_histogram;
		delete [] survival_histogram;

	}
}

AggregatedTemplateResult::AggregatedTemplateResult()
{
	image_number = -1;
	number_of_received_results = 0;
	total_number_of_ccs = 0.0f;

	collated_data_array = NULL;
	collated_mip_data = NULL;
	collated_psi_data = NULL;
	collated_theta_data = NULL;
	collated_phi_data = NULL;
	collated_defocus_data = NULL;
	collated_pixel_size_data = NULL;
	collated_pixel_sums = NULL;
	collated_pixel_square_sums = NULL;
	collated_histogram_data = NULL;
}

AggregatedTemplateResult::~AggregatedTemplateResult()
{
	if (collated_data_array != NULL) delete [] collated_data_array;

}

void AggregatedTemplateResult::AddResult(float *result_array, long array_size, int result_number, int number_of_expected_results)
{

	int offset = number_of_meta_data_values;

	if (collated_data_array == NULL)
	{
		collated_data_array = new float[array_size];
		ZeroFloatArray(collated_data_array, array_size);
		number_of_received_results = 0;
		total_number_of_ccs = 0.0f;

		// nasty..

		collated_mip_data 			= &collated_data_array[offset + int(result_array[2]) * 0];
		collated_psi_data 			= &collated_data_array[offset + int(result_array[2]) * 1];
		collated_theta_data 		= &collated_data_array[offset + int(result_array[2]) * 2];
		collated_phi_data 			= &collated_data_array[offset + int(result_array[2]) * 3];
		collated_defocus_data 		= &collated_data_array[offset + int(result_array[2]) * 4];
		collated_pixel_size_data	= &collated_data_array[offset + int(result_array[2]) * 5];
		collated_pixel_sums 		= &collated_data_array[offset + int(result_array[2]) * 6];
		collated_pixel_square_sums 	= &collated_data_array[offset + int(result_array[2]) * 7];

		collated_histogram_data = (long *) &collated_data_array[offset  + int(result_array[2]) * 8];

		collated_data_array[0] = result_array[0];
		collated_data_array[1] = result_array[1];
		collated_data_array[2] = result_array[2];
		collated_data_array[3] = result_array[3];

		collated_data_array[5] = result_array[5];

	}

	total_number_of_ccs += result_array[4];


	float *result_mip_data 			= &result_array[offset + int(result_array[2]) * 0];
	float *result_psi_data 			= &result_array[offset + int(result_array[2]) * 1];
	float *result_theta_data 		= &result_array[offset + int(result_array[2]) * 2];
	float *result_phi_data 			= &result_array[offset + int(result_array[2]) * 3];
	float *result_defocus_data 		= &result_array[offset + int(result_array[2]) * 4];
	float *result_pixel_size_data 	= &result_array[offset + int(result_array[2]) * 5];
	float *result_pixel_sums 		= &result_array[offset + int(result_array[2]) * 6];
	float *result_pixel_square_sums = &result_array[offset + int(result_array[2]) * 7];

	long *input_histogram_data = (long *) &result_array[offset + int(result_array[2]) * 8];

	long pixel_counter;
	long result_array_counter;

	// handle the images..

	for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
	{
		if (result_mip_data[pixel_counter] >  collated_mip_data[pixel_counter])
		{
			collated_mip_data[pixel_counter] = result_mip_data[pixel_counter];
			collated_psi_data[pixel_counter] = result_psi_data[pixel_counter];
			collated_theta_data[pixel_counter] = result_theta_data[pixel_counter];
			collated_phi_data[pixel_counter] = result_phi_data[pixel_counter];
			collated_defocus_data[pixel_counter] = result_defocus_data[pixel_counter];
			collated_pixel_size_data[pixel_counter] = result_pixel_size_data[pixel_counter];
		}
	}


	// sums and sum of squares

	for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
	{
		collated_pixel_sums[pixel_counter] += result_pixel_sums[pixel_counter];
	}

	for (pixel_counter = 0; pixel_counter <  int(result_array[2]); pixel_counter++)
	{
		collated_pixel_square_sums[pixel_counter] += result_pixel_square_sums[pixel_counter];
	}

	// handle the histogram..

	for (pixel_counter = 0; pixel_counter < histogram_number_of_points; pixel_counter++)
	{
		collated_histogram_data[pixel_counter] += input_histogram_data[pixel_counter];
	}

	number_of_received_results++;
	wxPrintf("Received %i of %i results\n", number_of_received_results, number_of_expected_results);
}
