#include "../../core/core_headers.h"


// Values for data that are passed around in the results.
const int number_of_output_images = 8; //mip, scaledmip, psi, theta, phi, pixel, defocus, sums, sqsums
const int number_of_meta_data_values = 6; // img_x, img_y, number cccs, histogram values.
const int MAX_ALLOWED_NUMBER_OF_PEAKS = 1000; // An error will be thrown and job aborted if this number of peaks is exceeded in the make template results block
const int blur_hack_n_frames = 128;
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

	bool do_shift_blur_hack = true;
	bool do_exposure_filter_hack = true;
	float* shift_hack_x;
	float* shift_hack_y;
	float* shift_hack_d;

	if (do_shift_blur_hack)
	{
		wxPrintf("Doing the shift blur hack\n");
		wxString tmp_string = input_search_images_filename;

		bool is_set = false;
		if (tmp_string.AfterLast('/').StartsWith("May06_12.48.59_20_1"))shift_hack_x = new float[blur_hack_n_frames]{0.179864,0.226035,0.221800,0.101352,-0.028329,-0.082596,-0.054399,-0.002695,0.053346,0.122100,0.190697,0.215507,0.225322,0.271974,0.290356,0.289413,0.314659,0.307740,0.315842,0.341318,0.332869,0.323045,0.315857,0.289356,0.283353,0.282720,0.294014,0.296343,0.289211,0.306148,0.309726,0.305844,0.323079,0.316819,0.293598,0.282043,0.260122,0.241315,0.248519,0.246461,0.228426,0.232809,0.257193,0.268968,0.248897,0.225143,0.219199,0.245235,0.249922,0.214677,0.192343,0.219558,0.237092,0.214007,0.175625,0.145734,0.149653,0.159645,0.162509,0.165749,0.202934,0.221146,0.203262,0.206492,0.210446,0.193382,0.159791,0.130018,0.126414,0.145725,0.132738,0.113221,0.119900,0.125603,0.115583,0.103640,0.101890,0.121894,0.132855,0.122629,0.101249,0.089128,0.123589,0.152442,0.150613,0.136710,0.131279,0.156506,0.163083,0.139541,0.116563,0.100642,0.108566,0.122189,0.117298,0.119576,0.127128,0.125971,0.135680,0.144747,0.115492,0.114149,0.104267,0.082058,0.106594,0.106722,0.097531,0.120948,0.112326,0.112273,0.117238,0.086524,0.075906,0.080209,0.060397,0.049039,0.072293,0.076256,0.042657,0.024812,0.044528,0.062999,0.054148,0.015800,0.043612,0.090049,0.084364,0.075686} ; shift_hack_y = new float[blur_hack_n_frames]{-0.598946,-0.618100,-0.596117,-0.617086,-0.672865,-0.625347,-0.538435,-0.430179,-0.334469,-0.265841,-0.194734,-0.131746,-0.059495,0.013713,0.037070,0.048825,0.076130,0.097069,0.117404,0.117561,0.101573,0.096835,0.093847,0.102045,0.111480,0.093396,0.085382,0.101272,0.109509,0.122307,0.107102,0.104257,0.136126,0.126461,0.108804,0.117392,0.103071,0.106639,0.115067,0.077577,0.087286,0.092084,0.073823,0.089815,0.085058,0.073099,0.088257,0.086413,0.083106,0.084351,0.080511,0.096522,0.104486,0.097786,0.090173,0.101890,0.120414,0.129854,0.125678,0.112223,0.105365,0.088581,0.074687,0.047516,0.022550,0.035412,0.037453,0.030237,0.053876,0.061945,0.052706,0.048347,0.036499,0.040646,0.063917,0.073599,0.062171,0.064935,0.056181,0.049095,0.064211,0.052632,0.043471,0.024568,0.005579,0.050303,0.035240,-0.015680,0.019802,0.035353,0.044421,0.053848,0.017128,0.016784,0.053898,0.035002,0.008670,0.011114,0.016531,0.046956,0.066306,0.051566,0.051177,0.073602,0.073497,0.084179,0.062365,0.039128,0.032026,0.006703,-0.009745,0.017017,0.008955,-0.021688,0.000307,-0.004085,0.014546,0.026009,0.017563,0.038441,0.074100,0.056276,0.012083,0.010376,0.022198,-0.028871,-0.047565,-0.005983} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.45.05_28_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.076964,-1.155506,-1.145443,-0.984338,-0.839038,-0.718226,-0.624298,-0.528425,-0.449828,-0.422612,-0.383440,-0.311139,-0.274371,-0.271291,-0.250447,-0.205356,-0.161654,-0.123947,-0.118827,-0.117010,-0.093291,-0.080379,-0.072606,-0.062191,-0.065524,-0.061106,-0.049512,-0.036031,-0.008826,0.000610,-0.011463,-0.013589,-0.005986,-0.001644,-0.005619,-0.028267,-0.042511,-0.020780,-0.001382,-0.020484,-0.040623,-0.040835,-0.019477,0.007584,0.011379,0.014358,0.036602,0.068958,0.093039,0.081050,0.037517,-0.001317,-0.021377,-0.010435,-0.017181,-0.023696,-0.005392,0.014390,0.041199,0.036931,0.020179,0.027502,0.028056,-0.000602,-0.016208,-0.006700,0.019262,0.026640,0.014726,0.018022,0.033737,0.049270,0.056434,0.050054,0.043404,0.061251,0.077409,0.073430,0.070401,0.069671,0.058883,0.056776,0.044296,0.022048,0.033525,0.044510,0.044878,0.059644,0.055503,0.041537,0.025324,0.005805,-0.003916,-0.011251,-0.016063,0.006118,0.031692,0.017215,0.012335,0.016600,-0.006035,-0.025158,-0.059940,-0.075247,-0.021282,-0.002001,-0.023564,-0.016369,0.002356,0.027234,0.036664,0.004202,-0.000324,0.032056,0.018923,0.008487,0.013584,0.009320,0.026066,0.030178,0.026082,0.032244,0.033297,0.018436,0.005931,-0.007085,-0.010846,-0.002416} ; shift_hack_y = new float[blur_hack_n_frames]{-1.467819,-1.517214,-1.509818,-1.429056,-1.370687,-1.277726,-1.175693,-1.080080,-0.991035,-0.904557,-0.822457,-0.751113,-0.695840,-0.643655,-0.589811,-0.561242,-0.522274,-0.480975,-0.470147,-0.447988,-0.423892,-0.404792,-0.369750,-0.361115,-0.370128,-0.344448,-0.323633,-0.328483,-0.331914,-0.337692,-0.335106,-0.309707,-0.291854,-0.314094,-0.314378,-0.290137,-0.288601,-0.285002,-0.282582,-0.281342,-0.260333,-0.243799,-0.251907,-0.256573,-0.241887,-0.230431,-0.234336,-0.220039,-0.208872,-0.199933,-0.172255,-0.177624,-0.196024,-0.200350,-0.201365,-0.189013,-0.185747,-0.190471,-0.184814,-0.171504,-0.132875,-0.125272,-0.139523,-0.125905,-0.127590,-0.122627,-0.113371,-0.149008,-0.158718,-0.150756,-0.160249,-0.136895,-0.136610,-0.151563,-0.127714,-0.135279,-0.130015,-0.112943,-0.138071,-0.139343,-0.120809,-0.122083,-0.110865,-0.116679,-0.125791,-0.107742,-0.104901,-0.117306,-0.136250,-0.125651,-0.118259,-0.108912,-0.108801,-0.124708,-0.117431,-0.107948,-0.123682,-0.113524,-0.107669,-0.132459,-0.134409,-0.108716,-0.077407,-0.065211,-0.055444,-0.054139,-0.046144,-0.047486,-0.075797,-0.077103,-0.061878,-0.080389,-0.108959,-0.111044,-0.072963,-0.059115,-0.082819,-0.100067,-0.098696,-0.080457,-0.061060,-0.081300,-0.097116,-0.085471,-0.065031,-0.051523,-0.060991,-0.060364} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.50.26_29_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.860774,-1.962393,-1.958552,-1.754390,-1.585977,-1.428430,-1.260191,-1.087192,-0.923631,-0.793156,-0.710755,-0.632715,-0.553784,-0.502921,-0.452333,-0.412102,-0.385443,-0.331616,-0.275966,-0.241083,-0.205945,-0.198026,-0.192671,-0.153718,-0.153161,-0.165344,-0.158132,-0.172187,-0.165460,-0.162236,-0.183776,-0.167540,-0.131628,-0.129621,-0.126603,-0.123831,-0.119153,-0.106729,-0.106861,-0.119296,-0.127083,-0.112719,-0.090237,-0.075550,-0.064035,-0.084219,-0.106067,-0.108499,-0.130085,-0.148391,-0.158799,-0.181190,-0.169487,-0.140168,-0.116247,-0.096757,-0.100750,-0.105727,-0.093690,-0.100167,-0.117723,-0.118864,-0.097247,-0.081332,-0.076198,-0.071591,-0.064151,-0.069055,-0.095217,-0.091904,-0.073798,-0.063896,-0.068488,-0.078849,-0.061827,-0.035449,-0.052047,-0.072663,-0.089552,-0.095539,-0.066959,-0.080564,-0.100587,-0.096776,-0.085914,-0.058776,-0.058550,-0.077307,-0.076176,-0.056590,-0.052555,-0.059115,-0.062732,-0.052961,-0.035063,-0.016675,-0.026943,-0.026780,-0.031456,-0.042350,-0.025431,-0.049961,-0.065025,-0.057886,-0.062645,-0.066251,-0.056546,-0.068219,-0.055251,-0.027388,-0.042629,-0.048347,-0.034883,-0.035197,-0.016551,-0.026348,-0.059994,-0.025986,-0.026078,-0.019345,-0.017727,-0.060651,-0.039760,-0.005751,-0.016681,-0.008185,0.001098,-0.007429} ; shift_hack_y = new float[blur_hack_n_frames]{-2.362815,-2.405574,-2.397419,-2.369353,-2.318621,-2.192277,-2.030024,-1.825586,-1.603813,-1.397378,-1.225573,-1.087353,-0.961576,-0.832365,-0.741603,-0.676229,-0.612156,-0.550275,-0.491437,-0.445003,-0.399288,-0.359874,-0.325914,-0.291133,-0.286091,-0.275162,-0.247801,-0.245983,-0.232563,-0.215429,-0.208775,-0.182715,-0.179121,-0.173521,-0.161846,-0.169603,-0.166114,-0.160631,-0.162945,-0.155614,-0.152739,-0.134762,-0.121702,-0.119578,-0.113578,-0.107618,-0.095590,-0.089334,-0.084852,-0.080421,-0.095136,-0.099723,-0.097732,-0.102186,-0.111258,-0.130107,-0.115832,-0.092733,-0.087592,-0.091352,-0.080471,-0.065658,-0.045364,-0.036452,-0.036104,-0.015140,-0.001442,-0.013132,-0.015564,-0.020044,-0.036823,-0.051894,-0.050212,-0.044861,-0.057083,-0.050912,-0.051354,-0.061376,-0.053407,-0.056939,-0.050603,-0.023510,-0.021304,-0.006618,0.019663,0.026852,0.007659,0.004076,-0.002199,-0.013679,-0.024244,-0.038140,-0.040906,-0.040361,-0.064303,-0.075792,-0.056251,-0.049573,-0.063932,-0.060643,-0.052560,-0.047403,-0.053413,-0.079241,-0.082263,-0.068710,-0.081464,-0.073058,-0.044059,-0.038331,-0.033911,-0.039224,-0.036667,-0.030730,-0.032876,-0.024582,-0.015947,-0.029322,-0.028846,-0.014901,0.006390,-0.013818,-0.045805,-0.006810,0.023244,0.024175,0.027946,0.020103} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.34.54_26_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.925367,-1.023360,-1.029377,-0.789848,-0.553345,-0.422511,-0.338979,-0.295183,-0.284493,-0.263596,-0.259635,-0.243867,-0.243310,-0.243398,-0.217928,-0.220929,-0.234173,-0.240186,-0.247569,-0.247493,-0.251521,-0.245186,-0.233937,-0.215851,-0.174106,-0.172916,-0.158652,-0.104544,-0.075578,-0.073348,-0.068141,-0.074826,-0.067497,-0.062865,-0.103411,-0.136533,-0.102523,-0.073498,-0.068008,-0.045348,-0.045347,-0.044355,-0.022644,-0.047190,-0.080025,-0.079503,-0.082949,-0.081670,-0.068697,-0.065058,-0.065765,-0.051288,-0.060814,-0.073237,-0.054142,-0.053357,-0.085401,-0.086186,-0.071696,-0.070841,-0.065233,-0.073974,-0.058881,-0.019729,-0.014816,-0.021567,-0.018349,-0.024217,-0.032488,-0.047774,-0.044070,-0.041435,-0.059851,-0.058041,-0.069456,-0.070796,-0.074941,-0.093839,-0.059830,-0.033890,-0.038582,-0.013069,0.005166,0.007310,0.003242,-0.024259,-0.040409,-0.025248,-0.035485,-0.028395,0.001982,-0.004035,-0.015172,-0.023542,-0.033759,-0.034003,-0.049783,-0.047660,0.002326,0.017911,0.032835,0.053116,0.062648,0.078402,0.066624,0.028526,0.007075,0.028061,0.033155,0.034264,0.044282,0.037257,0.032798,0.044916,0.021985,0.014987,0.023230,-0.009523,0.000005,0.029613,0.049949,0.051719,0.018207,0.017677,0.032202,-0.006135,-0.017790,0.003630} ; shift_hack_y = new float[blur_hack_n_frames]{-0.944960,-1.047369,-1.048236,-0.821022,-0.632156,-0.496834,-0.390889,-0.333466,-0.288562,-0.238561,-0.199403,-0.185534,-0.160070,-0.129960,-0.094711,-0.064952,-0.063925,-0.064606,-0.060546,-0.049025,-0.039550,-0.028910,-0.027580,-0.030181,-0.024559,-0.000989,-0.004457,-0.011646,0.001177,-0.005297,0.002212,0.000771,-0.005606,0.005308,0.007187,-0.004148,-0.023846,-0.022663,0.013564,0.007856,-0.024387,-0.028157,-0.022439,-0.000499,0.001265,-0.018326,-0.010222,0.005469,-0.007994,-0.023206,-0.017260,-0.014476,-0.011258,-0.008233,0.009394,0.036660,0.030886,0.020468,0.011494,-0.015163,-0.020307,-0.032419,-0.057815,-0.030260,-0.023854,-0.036643,-0.011780,-0.011910,0.000182,0.038539,0.030962,0.030429,0.057000,0.051648,0.064824,0.096843,0.086225,0.076743,0.095757,0.080328,0.053715,0.032226,-0.001772,0.005235,0.009561,-0.013301,0.002609,0.039765,0.060610,0.065763,0.045154,0.035735,0.062679,0.055709,0.052477,0.047079,0.022913,0.017335,0.024738,0.034429,0.048213,0.019714,0.002216,0.023389,0.043465,0.057155,0.013315,-0.027135,0.010552,0.081343,0.113979,0.094843,0.064254,0.091672,0.133195,0.090720,0.020444,0.011264,0.045169,0.044453,-0.001670,-0.014603,0.004882,0.046874,0.037788,0.024666,0.045475,0.042521} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_15.16.42_47_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.052099,-1.104188,-1.106793,-0.994128,-0.875416,-0.776361,-0.688781,-0.629467,-0.567376,-0.471045,-0.408957,-0.362317,-0.319491,-0.294355,-0.266817,-0.232852,-0.207430,-0.194119,-0.168796,-0.151539,-0.143000,-0.119391,-0.128435,-0.157505,-0.153080,-0.140019,-0.119477,-0.107872,-0.103830,-0.074583,-0.057543,-0.052851,-0.052967,-0.067161,-0.073717,-0.065454,-0.064100,-0.050578,-0.043938,-0.061586,-0.065787,-0.065468,-0.056630,-0.056490,-0.076281,-0.071640,-0.052799,-0.035343,-0.011483,-0.038468,-0.058225,-0.048609,-0.045636,-0.054137,-0.064528,-0.070168,-0.050251,-0.031754,-0.034524,-0.024595,-0.004387,0.008977,0.005297,0.017223,0.023533,0.002519,-0.007947,-0.003624,-0.008599,-0.011949,-0.012115,-0.025762,-0.009712,0.013489,-0.005848,-0.012895,-0.004751,-0.021517,-0.031495,-0.027485,-0.058240,-0.065003,-0.045241,-0.018580,-0.015609,-0.017928,-0.003737,-0.007714,-0.004451,-0.024393,-0.056922,-0.044304,-0.018821,-0.036510,-0.058394,-0.074923,-0.068132,-0.052247,-0.062506,-0.070201,-0.040394,-0.007119,0.017714,0.020886,0.008142,0.026401,0.033428,0.014045,0.005834,-0.024598,-0.028616,0.015443,0.013218,-0.008439,-0.009735,0.005414,0.025921,0.014068,-0.021782,-0.018851,0.004664,0.022810,0.026164,0.043329,0.063111,0.059465,0.056193,0.055328} ; shift_hack_y = new float[blur_hack_n_frames]{-0.851525,-0.911112,-0.900969,-0.798930,-0.706884,-0.603950,-0.538946,-0.473073,-0.371599,-0.292910,-0.249475,-0.206522,-0.163688,-0.133852,-0.112234,-0.091886,-0.091737,-0.074293,-0.064798,-0.061245,-0.025394,-0.001852,0.004550,0.027976,0.021355,0.000978,0.005850,-0.001426,-0.011208,-0.023253,-0.042684,-0.035370,-0.014224,-0.029532,-0.044951,-0.027068,-0.015952,-0.011837,-0.006072,-0.013504,-0.002449,0.007042,-0.011589,-0.004979,-0.010410,-0.020558,0.005220,-0.002825,-0.005116,0.019146,0.015376,0.019449,0.019553,-0.008157,-0.012305,-0.015200,-0.026907,-0.021098,-0.024895,-0.034552,-0.018035,0.000363,0.001768,-0.016115,-0.031780,-0.020942,-0.008344,-0.023129,-0.024545,-0.003909,-0.018875,-0.035585,-0.037265,-0.015029,0.003255,-0.017467,-0.031572,-0.010689,0.011463,0.018906,0.006295,-0.020328,-0.022160,-0.016350,-0.008311,-0.019735,-0.027144,-0.010619,-0.012796,-0.037074,-0.045404,-0.039025,-0.012765,-0.015931,-0.036010,0.008367,0.039230,0.022034,0.025854,0.002522,-0.012934,-0.021688,-0.054468,-0.031013,0.001560,-0.029359,-0.026102,-0.000642,-0.016200,-0.033973,-0.029919,-0.040498,-0.042071,-0.044720,-0.077861,-0.065827,-0.021644,-0.032628,-0.036057,-0.028376,-0.028550,0.012419,-0.016102,-0.068605,-0.052374,-0.061945,-0.081528,-0.057510} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_18.59.55_91_1"))shift_hack_x = new float[blur_hack_n_frames]{0.335093,0.360048,0.366411,0.284682,0.208394,0.190741,0.197390,0.220488,0.235825,0.252992,0.290978,0.311947,0.327797,0.332588,0.332479,0.340972,0.334571,0.309322,0.289002,0.288453,0.274272,0.262773,0.237156,0.222483,0.227828,0.217993,0.199316,0.194724,0.194123,0.206117,0.219078,0.229439,0.235757,0.217369,0.216893,0.212926,0.193421,0.176284,0.154940,0.153922,0.160943,0.155043,0.138120,0.133777,0.132472,0.112589,0.102600,0.098796,0.081896,0.080267,0.089023,0.096501,0.108804,0.114018,0.125042,0.134081,0.147436,0.155262,0.159842,0.175957,0.173353,0.155282,0.151481,0.157438,0.156293,0.136148,0.121574,0.124239,0.139187,0.157332,0.151528,0.146093,0.152878,0.141681,0.141431,0.147520,0.118465,0.088951,0.074765,0.094770,0.109607,0.099075,0.091665,0.103575,0.122996,0.132850,0.112386,0.089683,0.073774,0.075687,0.085538,0.062030,0.041398,0.039237,0.058018,0.069788,0.058693,0.038924,0.057885,0.060357,0.044461,0.040035,0.046922,0.054722,0.055828,0.044260,0.043697,0.051838,0.048610,0.047573,0.048571,0.063683,0.077029,0.074528,0.075370,0.089672,0.081183,0.054560,0.045074,0.051858,0.070236,0.075838,0.055319,0.058198,0.072676,0.064573,0.058204,0.061795} ; shift_hack_y = new float[blur_hack_n_frames]{-0.622593,-0.604605,-0.587701,-0.716842,-0.838798,-0.797828,-0.662467,-0.500800,-0.364136,-0.202665,-0.066392,0.009897,0.064545,0.114971,0.159751,0.204639,0.223119,0.232300,0.256529,0.285616,0.285728,0.275097,0.276256,0.254215,0.235761,0.237031,0.231301,0.228926,0.233882,0.228425,0.237159,0.231567,0.205703,0.203736,0.204244,0.181518,0.156446,0.157795,0.170916,0.172400,0.163582,0.150901,0.147325,0.157373,0.138352,0.110806,0.110999,0.098351,0.071645,0.065782,0.066933,0.069055,0.085876,0.072648,0.053027,0.065741,0.061438,0.041545,0.027120,0.003312,0.001500,0.009583,0.002737,-0.003149,-0.007708,-0.007878,-0.003989,-0.006516,-0.022692,-0.032992,-0.019770,-0.002692,-0.011438,-0.029941,-0.026706,-0.007082,0.004438,0.006767,-0.007812,-0.003497,0.015589,0.001662,-0.000338,-0.009469,-0.040919,-0.038758,-0.045091,-0.067462,-0.075862,-0.088811,-0.070313,-0.039538,-0.034290,-0.022505,-0.016649,-0.016387,-0.003673,-0.015853,-0.035514,-0.033767,-0.040905,-0.035778,-0.023472,-0.044306,-0.036701,-0.017075,-0.035044,-0.029138,-0.026442,-0.032861,-0.015115,-0.021501,-0.035962,-0.025558,-0.019390,-0.024638,-0.022623,-0.007725,-0.008051,-0.018726,-0.025332,-0.027195,-0.011372,-0.010132,-0.043486,-0.051634,-0.038828,-0.041469} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_22.57.31_135_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.726345,-0.817400,-0.827650,-0.592464,-0.371410,-0.219880,-0.137395,-0.123083,-0.035964,0.011742,0.019388,0.062800,0.053864,0.057609,0.107990,0.116823,0.119373,0.143540,0.122102,0.137352,0.168638,0.141041,0.111969,0.115204,0.118926,0.105922,0.075467,0.053733,0.075469,0.065481,0.022264,-0.015956,-0.008829,-0.010989,-0.087616,-0.155410,-0.152974,-0.137751,-0.170691,-0.195398,-0.201594,-0.198194,-0.164020,-0.195483,-0.223812,-0.194127,-0.204971,-0.198835,-0.197780,-0.217313,-0.185827,-0.188606,-0.214686,-0.214619,-0.234431,-0.210187,-0.193519,-0.216066,-0.218227,-0.206242,-0.190538,-0.186030,-0.184307,-0.164987,-0.162131,-0.149689,-0.120774,-0.117003,-0.099075,-0.095282,-0.111781,-0.101589,-0.104688,-0.114332,-0.112578,-0.126664,-0.121185,-0.115909,-0.120329,-0.105278,-0.109980,-0.096152,-0.062981,-0.063958,-0.062579,-0.077583,-0.093454,-0.063952,-0.056246,-0.086038,-0.099080,-0.090579,-0.060359,-0.060516,-0.079200,-0.056697,-0.025540,-0.024626,-0.033259,-0.025629,0.004127,0.008174,-0.025975,-0.021680,0.003505,0.010195,0.006373,-0.016505,-0.019473,0.002966,-0.009346,-0.001197,0.015907,0.015044,0.026102,0.024360,0.003104,0.006265,0.003239,-0.013854,-0.006609,-0.007059,-0.003516,0.009068,-0.001528,-0.015322,-0.008986,-0.006923} ; shift_hack_y = new float[blur_hack_n_frames]{-1.403664,-1.553470,-1.574262,-1.188093,-0.766339,-0.526527,-0.402431,-0.355127,-0.260675,-0.175913,-0.166417,-0.171900,-0.194220,-0.200146,-0.141691,-0.094322,-0.093506,-0.054652,-0.033318,0.021073,0.081527,0.033530,-0.013745,0.001920,-0.015277,-0.014317,-0.030747,-0.072245,-0.053742,-0.077068,-0.118957,-0.137913,-0.133881,-0.141186,-0.228846,-0.269325,-0.251938,-0.279488,-0.294687,-0.318848,-0.356699,-0.320020,-0.274602,-0.322445,-0.342677,-0.330358,-0.345015,-0.339360,-0.322686,-0.327461,-0.333186,-0.319665,-0.318442,-0.315701,-0.312963,-0.314525,-0.304139,-0.276511,-0.280074,-0.280981,-0.262831,-0.252463,-0.239221,-0.247579,-0.262660,-0.234002,-0.207582,-0.211858,-0.221385,-0.233628,-0.222064,-0.189846,-0.188508,-0.190701,-0.187171,-0.185874,-0.163570,-0.156561,-0.161438,-0.153859,-0.149623,-0.136155,-0.127966,-0.120200,-0.105952,-0.105404,-0.112736,-0.110592,-0.110750,-0.108453,-0.105614,-0.117986,-0.128893,-0.129502,-0.134946,-0.133076,-0.118521,-0.129703,-0.135414,-0.127147,-0.105228,-0.088561,-0.103347,-0.100784,-0.085580,-0.069008,-0.044537,-0.054044,-0.070554,-0.062247,-0.072665,-0.080894,-0.076397,-0.088944,-0.095755,-0.071329,-0.046217,-0.038292,-0.038569,-0.058710,-0.059841,-0.054419,-0.070240,-0.066704,-0.065482,-0.055437,-0.047046,-0.057536} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_19.29.44_96_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.667405,-0.775482,-0.767187,-0.534745,-0.319124,-0.179742,-0.096038,-0.037986,0.021073,0.038239,0.029811,0.058417,0.064018,0.058926,0.084903,0.074887,0.070877,0.100696,0.084311,0.092148,0.114372,0.113263,0.111034,0.121374,0.118709,0.107642,0.087566,0.083473,0.102889,0.114671,0.105530,0.097021,0.091601,0.075670,0.078291,0.081250,0.054157,0.049573,0.065429,0.055789,0.079262,0.104325,0.086619,0.082216,0.083113,0.079628,0.109268,0.123697,0.112803,0.129963,0.148129,0.118812,0.107034,0.123410,0.110817,0.107371,0.100155,0.077868,0.085259,0.079799,0.055377,0.049039,0.038371,0.037731,0.049401,0.076809,0.088477,0.078076,0.079039,0.083252,0.070550,0.050797,0.041950,0.044074,0.037283,0.025438,0.026722,0.015602,0.010111,0.019294,0.023282,0.030226,0.045650,0.040125,0.033493,0.046596,0.058810,0.051489,0.034094,0.017716,0.006803,0.018820,0.045456,0.036794,0.040355,0.064864,0.072911,0.064755,0.055158,0.062065,0.051151,0.032602,0.011162,0.013628,0.043347,0.058205,0.015371,-0.033738,-0.013446,0.030450,0.005642,-0.021406,-0.037352,-0.017881,0.032891,0.004534,-0.005078,0.030676,0.031726,0.032133,0.023866,-0.007967,-0.015314,-0.020761,-0.021398,-0.012226,-0.006943,-0.008176} ; shift_hack_y = new float[blur_hack_n_frames]{-0.736203,-0.882245,-0.873523,-0.546969,-0.271734,-0.081210,0.044140,0.104290,0.145846,0.147092,0.138860,0.159283,0.139472,0.121938,0.125230,0.108713,0.128651,0.136330,0.098322,0.093155,0.067704,0.065948,0.110936,0.102088,0.077039,0.071664,0.069981,0.105834,0.107241,0.071335,0.062435,0.071007,0.080476,0.082828,0.064348,0.055093,0.061786,0.043826,0.057181,0.050714,0.025324,0.050422,0.056140,0.020492,0.006206,-0.002589,-0.005492,0.004007,-0.010521,-0.007491,0.026286,0.032410,0.019299,0.028522,0.012375,0.000704,-0.014164,-0.046304,-0.036151,-0.028316,-0.043522,-0.054737,-0.062970,-0.046805,-0.032230,-0.039886,-0.049253,-0.036482,-0.004117,-0.009485,-0.043998,-0.034575,-0.029878,-0.040184,-0.036260,-0.055664,-0.054929,-0.044086,-0.068073,-0.072283,-0.054125,-0.073611,-0.065444,-0.046406,-0.047920,-0.053236,-0.061523,-0.079971,-0.078349,-0.076255,-0.080035,-0.083513,-0.084246,-0.066233,-0.058702,-0.076225,-0.092591,-0.092037,-0.103925,-0.095650,-0.102391,-0.087520,-0.063656,-0.069573,-0.068922,-0.047309,-0.040694,-0.054299,-0.083946,-0.093019,-0.073790,-0.094316,-0.110555,-0.120100,-0.100858,-0.060573,-0.041882,-0.057938,-0.074150,-0.045831,-0.030673,-0.040423,-0.048295,-0.060366,-0.047728,-0.037569,-0.042736,-0.046370} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_19.35.09_97_1"))shift_hack_x = new float[blur_hack_n_frames]{0.282880,0.286952,0.276426,0.299817,0.304548,0.279127,0.265679,0.235783,0.221279,0.201801,0.191177,0.194052,0.171191,0.134810,0.135471,0.142683,0.113660,0.079392,0.051446,0.072315,0.104166,0.083731,0.058721,0.073964,0.095828,0.111622,0.101914,0.081077,0.094408,0.094867,0.098807,0.084475,0.058882,0.059642,0.060420,0.034698,0.023660,0.021991,0.020496,0.034993,0.028614,0.009485,0.010734,0.008664,-0.000379,0.019107,0.027999,0.024220,0.026893,0.036400,0.052985,0.041243,-0.000188,-0.008819,-0.019317,-0.020796,-0.014962,-0.038841,-0.022291,0.015506,0.015190,0.020320,0.030451,0.026929,0.037925,0.027663,0.014615,0.027981,0.029341,0.021787,0.019311,0.014682,0.017069,0.013651,-0.006230,-0.003558,0.008193,0.019032,0.034311,0.044052,0.041773,0.027758,0.033912,0.042366,0.047113,0.050113,0.031365,0.035795,0.067445,0.046633,0.038444,0.040932,0.018323,0.035004,0.047979,0.038000,0.055513,0.046945,0.038595,0.056289,0.044926,0.023199,0.010574,0.024144,0.041824,0.037177,0.035186,0.043613,0.053812,0.054361,0.035718,0.034570,0.030369,0.021192,0.033201,0.051639,0.055707,0.041539,0.046384,0.045487,0.028711,0.044528,0.057376,0.076545,0.080666,0.081846,0.089563,0.076385} ; shift_hack_y = new float[blur_hack_n_frames]{0.208640,0.191309,0.200263,0.224292,0.227324,0.251435,0.284501,0.281041,0.302737,0.352398,0.399852,0.425843,0.432229,0.435573,0.450858,0.461590,0.440708,0.419098,0.398222,0.393627,0.408980,0.407247,0.393234,0.397549,0.378226,0.375515,0.369032,0.358562,0.365740,0.352843,0.319813,0.313680,0.306376,0.290370,0.282833,0.262121,0.259940,0.261719,0.247985,0.257384,0.263186,0.243683,0.229540,0.225254,0.225500,0.230358,0.218540,0.208234,0.192323,0.181736,0.184931,0.161210,0.142964,0.121733,0.094905,0.112671,0.129855,0.112574,0.110837,0.099646,0.099073,0.111367,0.094975,0.077630,0.078816,0.084780,0.104578,0.093852,0.075718,0.082507,0.081718,0.065707,0.049909,0.018498,0.032305,0.071895,0.051488,0.037911,0.057360,0.072330,0.075629,0.060784,0.036775,0.044325,0.047026,0.020460,0.001404,0.007915,0.036296,0.039365,0.015139,0.006756,0.018636,0.045357,0.060249,0.023195,-0.012598,0.011678,0.035693,0.039788,0.014436,-0.026249,-0.006950,0.016412,0.002962,-0.015743,-0.029788,-0.029857,-0.010451,-0.012969,-0.016718,-0.006955,0.010025,0.008554,0.009384,-0.006317,-0.047598,-0.019964,0.012216,-0.024595,-0.009852,-0.000087,0.012211,0.046409,-0.001666,-0.024317,-0.018141,-0.009417} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_20.34.46_108_1"))shift_hack_x = new float[blur_hack_n_frames]{1.358260,1.405086,1.399419,1.321078,1.248741,1.143784,1.049547,0.968444,0.875247,0.772522,0.688884,0.608526,0.532621,0.482482,0.454740,0.444263,0.426051,0.387376,0.370328,0.365753,0.342677,0.326084,0.325327,0.310464,0.294362,0.284012,0.276256,0.268254,0.236279,0.223606,0.214865,0.198512,0.196120,0.187294,0.174897,0.179711,0.151602,0.137608,0.163196,0.161847,0.152250,0.144606,0.144252,0.156396,0.152925,0.141814,0.141462,0.137212,0.126662,0.117160,0.110128,0.112115,0.107432,0.097239,0.089843,0.092037,0.087521,0.079041,0.074193,0.067565,0.072065,0.074454,0.060612,0.061257,0.069626,0.064114,0.052538,0.050758,0.062143,0.078610,0.087205,0.084122,0.079086,0.074239,0.065248,0.067422,0.066839,0.049917,0.047242,0.058919,0.074384,0.071366,0.058257,0.055240,0.063512,0.069079,0.067475,0.063777,0.059780,0.053833,0.060740,0.064800,0.049898,0.042732,0.054835,0.064454,0.063029,0.053977,0.059953,0.072800,0.068137,0.051299,0.056952,0.066265,0.057433,0.058338,0.056358,0.059413,0.071767,0.062463,0.063342,0.079595,0.074744,0.072433,0.067675,0.066570,0.075646,0.073895,0.063999,0.057658,0.050665,0.043063,0.031201,0.031421,0.039279,0.041971,0.045609,0.044178} ; shift_hack_y = new float[blur_hack_n_frames]{1.305334,1.330049,1.318686,1.307926,1.322193,1.284659,1.186285,1.061071,0.966502,0.884427,0.817917,0.759776,0.686114,0.639974,0.599183,0.561831,0.546582,0.505101,0.455352,0.420061,0.404221,0.403577,0.393753,0.368948,0.346134,0.335299,0.326989,0.304153,0.285684,0.251745,0.215217,0.214789,0.214229,0.203594,0.195267,0.183520,0.192156,0.195911,0.182829,0.180400,0.166828,0.150569,0.130814,0.101641,0.104693,0.108233,0.095714,0.104957,0.105836,0.110473,0.131534,0.113119,0.091962,0.088979,0.073921,0.067917,0.074746,0.064527,0.063445,0.072785,0.067999,0.067009,0.070272,0.068337,0.062083,0.063628,0.064062,0.059207,0.063855,0.071661,0.070141,0.059172,0.052332,0.057627,0.058712,0.066225,0.053730,0.024434,0.039688,0.056878,0.042619,0.039755,0.023440,0.022387,0.052532,0.052560,0.040541,0.050918,0.055624,0.066019,0.064958,0.043952,0.051793,0.054882,0.053177,0.068121,0.065946,0.080586,0.097924,0.086703,0.091888,0.092455,0.068722,0.059411,0.066017,0.069611,0.060757,0.054210,0.058742,0.060948,0.063993,0.054860,0.041492,0.041355,0.051920,0.040100,0.045681,0.060658,0.046193,0.030146,0.032989,0.035467,0.047549,0.065879,0.052205,0.077765,0.091742,0.071663} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_22.46.52_133_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.172645,-1.342608,-1.331696,-0.991606,-0.699194,-0.425268,-0.245793,-0.093280,0.040335,0.098677,0.148318,0.183739,0.204032,0.221149,0.203013,0.188542,0.214307,0.216180,0.197999,0.167308,0.157028,0.183773,0.186421,0.172240,0.185256,0.186262,0.185001,0.169554,0.164032,0.172675,0.166181,0.141329,0.121173,0.130064,0.154644,0.159859,0.176188,0.200665,0.191580,0.197048,0.218833,0.217226,0.186028,0.150831,0.120449,0.134813,0.135453,0.091232,0.094753,0.117514,0.110890,0.137573,0.170117,0.159636,0.184383,0.181551,0.167094,0.186168,0.158630,0.133877,0.157904,0.156737,0.162048,0.146995,0.113321,0.133788,0.121269,0.084357,0.087356,0.088318,0.093266,0.092065,0.053476,0.038026,0.041979,0.032026,0.010443,0.043411,0.082016,0.076663,0.081730,0.064760,0.037010,0.046651,0.010344,-0.011785,0.018605,0.034419,0.064059,0.051852,0.050221,0.074973,0.045397,0.034952,0.066242,0.061980,0.040323,0.022420,0.004755,0.020824,0.030591,-0.015291,-0.007221,0.051846,0.049655,0.004310,-0.016515,-0.001215,0.016570,-0.028614,-0.050539,-0.009357,0.040569,0.030734,-0.010272,-0.010870,0.035848,0.028933,-0.014560,-0.026522,-0.009850,0.006337,-0.001096,-0.019345,-0.016885,0.016273,-0.001267,-0.002901} ; shift_hack_y = new float[blur_hack_n_frames]{-0.949452,-1.067334,-1.054373,-0.816981,-0.572542,-0.381067,-0.283475,-0.204021,-0.135737,-0.078275,-0.037338,-0.027804,-0.017854,0.011109,0.015388,0.009809,-0.002214,-0.010735,0.016368,0.012873,0.004390,0.033485,0.063150,0.082577,0.073974,0.053951,0.055837,0.064823,0.051916,0.022926,0.008230,0.018948,0.002249,-0.012682,-0.005349,0.027855,0.023254,-0.007900,-0.018789,-0.006244,0.013670,0.000750,-0.016484,0.009745,0.023559,0.009565,0.005405,-0.014161,-0.014829,-0.031526,-0.048297,-0.037131,-0.015471,-0.031923,-0.036650,0.001057,0.008398,-0.022805,-0.027794,-0.012919,-0.003792,0.018394,-0.014471,-0.038928,-0.013114,-0.015348,-0.045901,-0.049575,-0.063630,-0.086650,-0.067539,-0.051378,-0.045569,-0.039769,-0.048074,-0.057000,-0.038184,-0.031097,-0.035611,-0.058356,-0.062939,-0.051752,-0.069338,-0.049210,-0.050010,-0.093332,-0.091670,-0.090649,-0.088313,-0.052216,-0.081726,-0.110825,-0.087189,-0.075373,-0.056658,-0.071069,-0.108453,-0.101039,-0.088236,-0.112253,-0.125752,-0.141818,-0.091288,-0.042441,-0.086910,-0.097821,-0.065562,-0.047123,-0.060120,-0.094844,-0.106683,-0.074544,-0.058432,-0.051510,-0.069620,-0.049231,-0.051862,-0.095952,-0.070800,-0.074396,-0.073017,-0.024497,-0.033525,-0.009562,0.005416,-0.035664,-0.029579,-0.026852} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_16.17.41_59_1"))shift_hack_x = new float[blur_hack_n_frames]{-2.329161,-2.510585,-2.479475,-2.181293,-1.950779,-1.668302,-1.394380,-1.139286,-0.916945,-0.737014,-0.618840,-0.552155,-0.468922,-0.369114,-0.306351,-0.265349,-0.245306,-0.215485,-0.142964,-0.109324,-0.125118,-0.119350,-0.101063,-0.060571,-0.041776,-0.069729,-0.086292,-0.077120,-0.065230,-0.071450,-0.079295,-0.070887,-0.061287,-0.037153,-0.015752,-0.009026,0.000665,-0.011395,-0.031049,-0.024330,-0.016100,-0.021838,-0.029966,-0.033762,-0.029569,-0.014376,-0.009310,-0.026765,-0.015595,-0.011320,-0.007365,0.031451,0.047299,0.057634,0.063173,0.047888,0.065781,0.067123,0.039644,0.027443,0.011756,0.021399,0.034955,0.018568,0.013273,0.016681,0.017834,0.020199,0.030329,0.024335,0.026316,0.039508,0.040307,0.037149,0.030353,0.034128,0.047819,0.033937,0.014273,0.005392,0.011648,0.040067,0.025142,0.006545,0.019259,0.021089,0.034156,0.047453,0.012978,0.014628,0.031296,0.019959,0.032325,0.025020,0.021922,0.048127,0.044712,0.044592,0.051924,0.043634,0.053342,0.050034,0.020239,0.007182,-0.000498,-0.002518,0.024101,0.026100,0.021092,0.051800,0.060641,0.062827,0.060563,0.045352,0.049649,0.038086,0.042303,0.032320,0.015357,0.046444,0.056651,0.055299,0.063693,0.038058,0.040125,0.042293,0.022709,0.031107} ; shift_hack_y = new float[blur_hack_n_frames]{-6.238672,-6.631771,-6.540968,-5.956625,-5.546288,-4.916944,-4.145427,-3.392611,-2.851956,-2.367688,-1.975830,-1.654014,-1.359872,-1.160797,-0.979304,-0.829865,-0.695774,-0.555724,-0.454846,-0.373588,-0.307991,-0.249691,-0.193358,-0.175031,-0.158955,-0.129011,-0.100249,-0.055491,-0.044981,-0.025633,0.007832,0.005220,-0.007799,-0.002258,0.002514,0.001092,0.005298,0.022701,0.056633,0.089432,0.085794,0.069040,0.099860,0.117727,0.100071,0.094023,0.093721,0.100852,0.111989,0.104802,0.101039,0.095645,0.093565,0.106812,0.104707,0.093559,0.098871,0.102837,0.103501,0.115559,0.112015,0.110324,0.117097,0.119209,0.123006,0.136853,0.124893,0.114085,0.131312,0.130394,0.119909,0.112736,0.086699,0.071066,0.064531,0.075767,0.094250,0.076564,0.079177,0.085678,0.097887,0.109860,0.085394,0.057759,0.036278,0.032688,0.060234,0.045672,0.040099,0.044124,0.037973,0.066500,0.059301,0.039376,0.032800,0.014045,0.009948,0.034634,0.048639,0.042750,0.044014,0.042323,0.046463,0.040774,0.046387,0.045338,0.038651,0.043760,0.039392,0.026928,0.027121,0.009239,0.006572,0.019291,0.005048,0.024330,0.035401,0.019194,0.006594,-0.003675,0.001513,0.021410,0.006830,-0.005737,0.020764,0.023857,0.009860,0.017664} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_20.08.22_103_1"))shift_hack_x = new float[blur_hack_n_frames]{-2.653414,-3.371696,-3.016379,-1.983428,-1.923842,-1.694024,-1.426651,-1.077187,-0.846899,-2.109006,-3.257958,-3.044918,-2.926828,-2.836345,-1.431031,-1.361190,-1.251878,0.145692,0.116013,0.130698,0.194653,0.196211,0.186323,0.185209,0.174496,0.228351,0.254088,0.220209,0.239491,0.286109,0.335207,0.341416,0.332274,0.348861,0.369384,0.367614,0.344879,0.328191,0.361692,0.384744,0.378542,0.376136,0.376005,0.407344,0.399503,0.357515,0.335084,0.324961,0.332343,0.360090,0.363758,0.360550,0.374306,0.398763,0.412594,0.409267,0.391927,0.366823,0.365454,0.359523,0.340482,0.335603,0.328728,0.309971,0.308018,0.308261,0.310430,0.316664,0.285666,0.241968,0.242842,0.258015,0.251928,0.249277,0.251689,0.271205,0.268223,0.242086,0.225028,0.202781,0.183962,0.172110,0.153759,0.155353,0.164746,0.153614,0.159000,0.147170,0.112943,0.095803,0.108179,0.110452,0.106598,0.102291,0.096094,0.106137,0.112071,0.106965,0.100537,0.101514,0.090262,0.069848,0.072892,0.086519,0.078551,0.074381,0.071634,0.060308,0.097468,0.112653,0.104353,0.109896,0.102101,0.066900,0.066202,0.072423,0.046214,0.044929,0.056004,0.074742,0.088282,0.068569,0.052140,0.074500,0.069281,0.058455,0.064435,0.065483} ; shift_hack_y = new float[blur_hack_n_frames]{-2.558412,-3.163280,-2.851393,-2.034921,-2.047948,-1.840265,-1.577352,-1.207987,-0.936356,-1.921108,-2.831281,-2.652021,-2.524229,-2.373098,-1.154101,-1.140204,-1.086831,0.127418,0.159883,0.187764,0.229430,0.280537,0.319299,0.298260,0.272341,0.295991,0.288971,0.226572,0.212735,0.265477,0.298785,0.287354,0.292860,0.299984,0.316920,0.324155,0.289222,0.276864,0.294174,0.285897,0.282111,0.302859,0.315564,0.330343,0.333018,0.316426,0.317015,0.321163,0.307929,0.293351,0.268566,0.283529,0.317027,0.301408,0.290967,0.277563,0.288391,0.321789,0.298773,0.261520,0.261564,0.263696,0.261083,0.247183,0.222219,0.199838,0.205481,0.230026,0.220220,0.222286,0.227708,0.207028,0.208265,0.195278,0.173072,0.172137,0.155176,0.149096,0.158008,0.165553,0.193047,0.194565,0.166594,0.171037,0.171926,0.160294,0.138718,0.115217,0.108589,0.108687,0.113625,0.119098,0.119156,0.120622,0.093386,0.076359,0.102310,0.096045,0.072587,0.060655,0.059861,0.088276,0.101509,0.073962,0.069720,0.063835,0.050092,0.060534,0.070993,0.042501,0.031626,0.039264,0.034369,0.052634,0.069710,0.065929,0.064105,0.051057,0.056561,0.089938,0.077343,0.054899,0.036139,0.041504,0.078232,0.056984,0.038428,0.054697} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_20.13.40_104_1"))shift_hack_x = new float[blur_hack_n_frames]{0.083284,0.106872,0.111237,0.049423,-0.068108,-0.142842,-0.116819,-0.084634,-0.074175,-0.062029,-0.046381,-0.000347,0.057248,0.093066,0.127520,0.147441,0.143283,0.124464,0.111705,0.125548,0.111131,0.070572,0.064314,0.076053,0.084173,0.093399,0.080705,0.070077,0.064149,0.054178,0.054565,0.053553,0.033680,0.011666,0.010627,0.019668,0.009195,-0.010324,-0.008355,0.002646,0.015045,0.013636,0.013635,0.021111,0.023005,0.016877,0.016635,0.024547,0.015302,0.005181,0.021498,0.052316,0.068050,0.063480,0.038801,0.041668,0.043829,0.026880,0.028826,0.029800,0.035249,0.028112,0.031374,0.061257,0.074839,0.060506,0.050861,0.049265,0.060413,0.047457,0.018037,0.024443,0.023563,0.019823,0.036722,0.041502,0.036859,0.032886,0.014191,0.013183,0.017680,0.013453,0.024065,0.037380,0.043299,0.047253,0.041248,0.024317,0.022771,0.024353,0.022334,0.031388,0.040239,0.049760,0.050236,0.039529,0.035940,0.020231,0.003214,0.004146,0.014884,0.038558,0.036987,0.031736,0.062816,0.089131,0.091256,0.064899,0.061482,0.071120,0.051988,0.042919,0.037008,0.021210,0.022554,0.021037,0.019758,0.030380,0.041348,0.047247,0.061783,0.064240,0.058198,0.068671,0.057117,0.036877,0.036633,0.041380} ; shift_hack_y = new float[blur_hack_n_frames]{-1.009315,-0.998596,-0.984757,-1.086564,-1.237015,-1.239719,-1.130506,-0.973290,-0.799199,-0.664130,-0.571686,-0.466022,-0.355565,-0.251804,-0.179693,-0.164102,-0.123745,-0.109037,-0.104544,-0.074674,-0.074214,-0.086799,-0.075023,-0.059937,-0.048447,-0.036831,-0.034689,-0.040636,-0.040414,-0.036760,-0.033647,-0.005567,0.004663,-0.027079,-0.049057,-0.035679,-0.030826,-0.044848,-0.060938,-0.081367,-0.075155,-0.043607,-0.047732,-0.061207,-0.062581,-0.071299,-0.069798,-0.061498,-0.089949,-0.104430,-0.079430,-0.095956,-0.125889,-0.114057,-0.103747,-0.090962,-0.096492,-0.124774,-0.103709,-0.079170,-0.099224,-0.105729,-0.107438,-0.113966,-0.105392,-0.108334,-0.117386,-0.118094,-0.118301,-0.112439,-0.091869,-0.084633,-0.094332,-0.090386,-0.072629,-0.071994,-0.078184,-0.084521,-0.088270,-0.071200,-0.068950,-0.079100,-0.069852,-0.074263,-0.082292,-0.082733,-0.084161,-0.067042,-0.053688,-0.050928,-0.047808,-0.047101,-0.045764,-0.049928,-0.054497,-0.056855,-0.066093,-0.062658,-0.053679,-0.053888,-0.042046,-0.038470,-0.043030,-0.039773,-0.024740,-0.004185,-0.011548,-0.037043,-0.052347,-0.039328,-0.044105,-0.066373,-0.077972,-0.052352,-0.026880,-0.017685,-0.026394,-0.048414,-0.041537,-0.027550,-0.048448,-0.061723,-0.059794,-0.054701,-0.024877,-0.023584,-0.030890,-0.030981} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_13.55.43_30_1"))shift_hack_x = new float[blur_hack_n_frames]{-2.128584,-2.259270,-2.255443,-1.996976,-1.766649,-1.524842,-1.331799,-1.188035,-1.052832,-0.912493,-0.772897,-0.690421,-0.649310,-0.579585,-0.501470,-0.469086,-0.424069,-0.352555,-0.318817,-0.314178,-0.288765,-0.261919,-0.228214,-0.227436,-0.229851,-0.210923,-0.164262,-0.142239,-0.151886,-0.130606,-0.131153,-0.135931,-0.115138,-0.114740,-0.110843,-0.076165,-0.052329,-0.035038,-0.036957,-0.041238,-0.032478,-0.033194,-0.038658,-0.068656,-0.061372,-0.023526,-0.061001,-0.076647,-0.035624,-0.046560,-0.050116,-0.051245,-0.080593,-0.064247,-0.046556,-0.078784,-0.072332,-0.069879,-0.105316,-0.083305,-0.066170,-0.097722,-0.097110,-0.084330,-0.082500,-0.066522,-0.091485,-0.087492,-0.043052,-0.043585,-0.059956,-0.048994,-0.065760,-0.067783,-0.052085,-0.073750,-0.085434,-0.091447,-0.095234,-0.054447,-0.040288,-0.074040,-0.069330,-0.076892,-0.092100,-0.065239,-0.061758,-0.084328,-0.094676,-0.090656,-0.050140,-0.029801,-0.074182,-0.105160,-0.083431,-0.064506,-0.074584,-0.106880,-0.118197,-0.079650,-0.074076,-0.096996,-0.085795,-0.073332,-0.052843,-0.049354,-0.073505,-0.045453,-0.000268,-0.014156,-0.047418,-0.053353,-0.038162,-0.034785,-0.044374,-0.071403,-0.067091,-0.036577,-0.052868,-0.084357,-0.053991,-0.010012,-0.015539,0.006348,0.009904,0.022402,0.041240,0.021382} ; shift_hack_y = new float[blur_hack_n_frames]{-2.394580,-2.480039,-2.445298,-2.343915,-2.325697,-2.185215,-1.959916,-1.764920,-1.604015,-1.425106,-1.240183,-1.071856,-0.910645,-0.806177,-0.728031,-0.635456,-0.605760,-0.571895,-0.516694,-0.502584,-0.448263,-0.365870,-0.328982,-0.275241,-0.226550,-0.229846,-0.217000,-0.220454,-0.236683,-0.221587,-0.206162,-0.210654,-0.181296,-0.165381,-0.173041,-0.164377,-0.162006,-0.170689,-0.174669,-0.170386,-0.140034,-0.108063,-0.120103,-0.116569,-0.112690,-0.104390,-0.084379,-0.084495,-0.099385,-0.080932,-0.062085,-0.043521,-0.029914,-0.042576,-0.041698,-0.017959,0.007625,0.022868,0.026691,0.023826,0.023458,0.019327,0.005548,-0.000084,0.009689,0.001173,0.015426,0.037303,0.026534,0.036860,0.042474,0.034928,0.038737,0.037146,0.056247,0.083157,0.085367,0.082763,0.062512,0.067437,0.067608,0.042504,0.059822,0.067300,0.047908,0.057674,0.078337,0.084678,0.082562,0.058685,0.059079,0.087350,0.087333,0.064606,0.068592,0.070463,0.057502,0.056611,0.037874,0.027260,0.032048,0.032169,0.051828,0.061757,0.046620,0.047768,0.065600,0.080944,0.069666,0.065656,0.071287,0.093196,0.111499,0.105351,0.125494,0.107379,0.094797,0.122405,0.093309,0.070236,0.086380,0.055590,0.075166,0.073533,0.049757,0.086056,0.082479,0.077328} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_22.23.11_130_1"))shift_hack_x = new float[blur_hack_n_frames]{0.531525,0.580900,0.570100,0.474237,0.445523,0.387184,0.324620,0.313295,0.270564,0.231598,0.230292,0.225128,0.223274,0.236373,0.212476,0.183915,0.186425,0.188811,0.159008,0.148563,0.122378,0.076512,0.081944,0.121732,0.119086,0.117967,0.109141,0.075757,0.095366,0.106095,0.076052,0.053818,0.042713,0.051935,0.081353,0.076095,0.074139,0.068617,0.067611,0.085233,0.079234,0.057625,0.046920,0.046103,0.049166,0.039029,0.016342,0.032132,0.058850,0.057345,0.048529,0.046060,0.083045,0.087214,0.039946,0.029266,0.033621,0.049004,0.055615,0.027255,0.020626,0.037953,0.046403,0.068638,0.049459,0.020424,0.032035,0.040001,0.043817,0.039665,0.042210,0.028158,0.020482,0.024306,0.032405,0.043252,0.020614,0.004120,0.017212,0.022263,0.020336,0.027613,0.006334,0.009200,0.027411,0.038925,0.072558,0.073903,0.036131,0.045127,0.067014,0.050485,0.057241,0.049633,0.057902,0.076568,0.080682,0.097798,0.105734,0.071488,0.037187,0.025337,0.042188,0.033694,0.027934,0.023040,0.008210,0.025285,0.012760,0.014881,0.027627,0.001682,0.004125,0.027162,0.017990,0.021986,0.011426,-0.000201,-0.023915,-0.009420,0.004270,0.009127,0.027379,0.012216,0.029638,0.061421,0.062427,0.048841} ; shift_hack_y = new float[blur_hack_n_frames]{-0.700358,-0.670385,-0.658037,-0.781746,-0.897382,-0.889902,-0.823684,-0.767489,-0.674294,-0.557869,-0.485707,-0.428138,-0.385565,-0.340036,-0.298965,-0.270593,-0.257583,-0.248877,-0.223533,-0.181667,-0.174255,-0.180557,-0.188060,-0.168631,-0.140680,-0.140087,-0.149793,-0.141409,-0.111016,-0.097309,-0.097467,-0.094209,-0.100748,-0.101159,-0.110603,-0.123103,-0.122699,-0.119229,-0.094610,-0.089538,-0.096267,-0.095135,-0.094497,-0.070466,-0.060076,-0.098434,-0.095228,-0.055778,-0.046677,-0.056412,-0.040469,-0.050671,-0.041898,-0.039217,-0.085993,-0.051416,-0.017254,-0.062852,-0.088562,-0.091880,-0.076579,-0.070870,-0.116616,-0.106046,-0.083915,-0.085608,-0.084078,-0.097419,-0.085583,-0.069518,-0.085623,-0.087598,-0.065968,-0.052057,-0.054517,-0.068624,-0.078875,-0.093171,-0.116914,-0.130220,-0.121170,-0.126455,-0.132807,-0.106223,-0.107169,-0.114003,-0.104978,-0.126322,-0.108165,-0.064689,-0.096480,-0.093398,-0.040432,-0.042228,-0.061066,-0.079091,-0.092848,-0.087349,-0.078547,-0.080821,-0.090208,-0.075526,-0.040202,-0.058056,-0.089174,-0.102124,-0.118511,-0.103393,-0.106166,-0.131095,-0.102013,-0.086380,-0.086414,-0.060973,-0.095363,-0.116815,-0.102636,-0.113142,-0.069690,-0.048567,-0.091353,-0.075501,-0.042156,-0.038632,-0.052424,-0.076150,-0.074511,-0.065346} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_23.12.42_138_1"))shift_hack_x = new float[blur_hack_n_frames]{-0.894239,-0.881506,-0.864724,-0.975693,-1.113738,-1.122720,-1.023114,-0.891643,-0.772309,-0.638840,-0.524054,-0.435592,-0.345563,-0.256818,-0.195706,-0.147018,-0.105098,-0.071504,-0.073259,-0.056623,-0.028229,-0.013456,-0.012503,-0.033131,-0.026938,0.023415,0.037512,0.021495,0.026402,0.031806,0.049032,0.054478,0.053478,0.052308,0.036279,0.030439,0.049182,0.048124,0.036679,0.016890,0.023143,0.067494,0.081567,0.073573,0.068669,0.061958,0.061542,0.054122,0.031997,0.002851,0.002165,0.027087,0.033516,0.028024,0.036813,0.054722,0.075620,0.069744,0.059008,0.057311,0.048819,0.051432,0.062210,0.065769,0.059396,0.059820,0.052829,0.067073,0.080340,0.056827,0.044572,0.055252,0.049606,0.058334,0.058348,0.035775,0.045387,0.065884,0.074227,0.078501,0.074369,0.063384,0.067013,0.069183,0.057006,0.046154,0.048198,0.060232,0.064693,0.068660,0.052171,0.031408,0.031728,0.030491,0.019324,0.021254,0.010056,0.010060,0.024777,0.012505,0.017098,0.031251,0.034840,0.026364,0.030018,0.061940,0.080291,0.062000,0.049204,0.041729,0.074835,0.082162,0.035884,0.024699,0.011320,0.014683,0.033062,-0.004451,-0.003169,0.027421,0.038211,0.059711,0.045193,0.038914,0.041258,0.019863,0.013175,0.021905} ; shift_hack_y = new float[blur_hack_n_frames]{-1.952593,-1.965846,-1.925356,-2.040171,-2.228294,-2.176815,-1.914235,-1.622159,-1.420263,-1.214512,-0.982982,-0.771413,-0.631448,-0.543916,-0.427556,-0.294835,-0.209375,-0.152004,-0.140127,-0.084797,-0.019460,-0.024034,-0.030591,-0.021870,-0.019026,-0.010321,-0.012340,-0.006628,0.017169,0.018943,0.018284,0.010983,0.034833,0.054872,0.046460,0.042243,0.028472,0.038640,0.060694,0.038329,0.038066,0.042139,0.031631,0.047476,0.040375,0.036527,0.046389,0.043559,0.071109,0.079764,0.064145,0.064204,0.052051,0.046578,0.040518,0.014917,0.006725,0.026656,0.043059,0.040634,0.055956,0.078972,0.074034,0.076813,0.083749,0.087269,0.081119,0.069849,0.055260,0.062176,0.094208,0.094768,0.072587,0.074628,0.076644,0.068847,0.053921,0.021661,0.015561,0.026272,0.021762,0.025295,0.030678,0.020200,0.023877,0.037563,0.033309,0.040904,0.029328,0.027855,0.063546,0.057045,0.036119,0.038555,0.044650,0.051502,0.036247,0.017997,0.014663,0.016609,0.022100,0.004988,0.006656,0.017539,0.015778,0.038811,0.036429,0.036516,0.043118,0.036330,0.047786,0.038366,0.019285,0.030847,0.022534,0.024516,0.043865,0.032095,0.041774,0.046735,0.031167,0.057286,0.094717,0.081142,0.065422,0.083463,0.086576,0.072998} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;
		if (tmp_string.AfterLast('/').StartsWith("May06_16.58.30_67_1"))shift_hack_x = new float[blur_hack_n_frames]{-1.170469,-1.210852,-1.196806,-1.154858,-1.132164,-1.032145,-0.897938,-0.810500,-0.729925,-0.609118,-0.518888,-0.435432,-0.375085,-0.347178,-0.291034,-0.239878,-0.206070,-0.170992,-0.162226,-0.144863,-0.133453,-0.130393,-0.102933,-0.111189,-0.119491,-0.104522,-0.107538,-0.082150,-0.060697,-0.082705,-0.096241,-0.085714,-0.072813,-0.065765,-0.077878,-0.107223,-0.104839,-0.082162,-0.089981,-0.117752,-0.124654,-0.101426,-0.072929,-0.074042,-0.080481,-0.058344,-0.037147,-0.029157,-0.029059,-0.028632,-0.019034,-0.014289,-0.017748,-0.006375,0.004973,-0.007154,-0.009563,0.005236,0.012431,0.010105,-0.003686,-0.003324,0.008172,0.016353,0.016215,0.012814,0.026179,0.040461,0.028581,0.025057,0.021967,-0.002459,-0.000467,0.006167,-0.010247,-0.010888,0.002732,0.014191,0.035161,0.033329,0.023581,0.035634,0.037321,0.023192,0.019034,0.014304,0.012130,0.016479,0.029157,0.034285,0.025228,0.028442,0.043970,0.055081,0.052500,0.031430,0.021355,0.014990,0.002667,-0.000640,-0.002812,-0.003238,0.011654,0.022309,0.029158,0.040728,0.020963,0.003460,0.003054,0.004842,0.012679,0.023786,0.028318,0.041629,0.039421,0.025500,0.029876,0.020620,0.013672,0.012691,0.000760,0.008972,0.023577,0.023383,0.036989,0.047363,0.043825,0.041521} ; shift_hack_y = new float[blur_hack_n_frames]{-2.724787,-2.766529,-2.736019,-2.758398,-2.781773,-2.616501,-2.420876,-2.193379,-1.921038,-1.671321,-1.475414,-1.298691,-1.148922,-1.028029,-0.910824,-0.816648,-0.773145,-0.730381,-0.664624,-0.606776,-0.549402,-0.514371,-0.522937,-0.492499,-0.452197,-0.464844,-0.472497,-0.450839,-0.444025,-0.444502,-0.448862,-0.458959,-0.434231,-0.423699,-0.426918,-0.419390,-0.406562,-0.380095,-0.347098,-0.337752,-0.338418,-0.340582,-0.335156,-0.306736,-0.283265,-0.282046,-0.290739,-0.274339,-0.283855,-0.315079,-0.315734,-0.313916,-0.320155,-0.297769,-0.297430,-0.298881,-0.266894,-0.256896,-0.262123,-0.254775,-0.257391,-0.256016,-0.229822,-0.213800,-0.209703,-0.207546,-0.212830,-0.211684,-0.208350,-0.199530,-0.181448,-0.184724,-0.177732,-0.155758,-0.156770,-0.137576,-0.121522,-0.149045,-0.159942,-0.156512,-0.147099,-0.135096,-0.144586,-0.168871,-0.169987,-0.149997,-0.144295,-0.153275,-0.153512,-0.131432,-0.111329,-0.103873,-0.108271,-0.101741,-0.095659,-0.084210,-0.090998,-0.106481,-0.103111,-0.099153,-0.108868,-0.105604,-0.110998,-0.111127,-0.100774,-0.098041,-0.075959,-0.073391,-0.072857,-0.080089,-0.093987,-0.088856,-0.088336,-0.103057,-0.077468,-0.066844,-0.072575,-0.056726,-0.037100,-0.018400,-0.035769,-0.045780,-0.038737,-0.046237,-0.068080,-0.099849,-0.098671,-0.086243} ; shift_hack_d = new float[blur_hack_n_frames]{1.155000,1.925000,2.695000,3.465000,4.235000,5.005000,5.775000,6.545000,7.315000,8.085000,8.855000,9.625000,10.395000,11.165000,11.935000,12.705000,13.475000,14.245000,15.015000,15.785000,16.555000,17.325000,18.095000,18.865000,19.635000,20.405000,21.175000,21.945000,22.715000,23.485000,24.255000,25.025000,25.795000,26.565000,27.335000,28.105000,28.875000,29.645000,30.415000,31.185000,31.955000,32.725000,33.495000,34.265000,35.035000,35.805000,36.575000,37.345000,38.115000,38.885000,39.655000,40.425000,41.195000,41.965000,42.735000,43.505000,44.275000,45.045000,45.815000,46.585000,47.355000,48.125000,48.895000,49.665000,50.435000,51.205000,51.975000,52.745000,53.515000,54.285000,55.055000,55.825000,56.595000,57.365000,58.135000,58.905000,59.675000,60.445000,61.215000,61.985000,62.755000,63.525000,64.295000,65.065000,65.835000,66.605000,67.375000,68.145000,68.915000,69.685000,70.455000,71.225000,71.995000,72.765000,73.535000,74.305000,75.075000,75.845000,76.615000,77.385000,78.155000,78.925000,79.695000,80.465000,81.235000,82.005000,82.775000,83.545000,84.315000,85.085000,85.855000,86.625000,87.395000,88.165000,88.935000,89.705000,90.475000,91.245000,92.015000,92.785000,93.555000,94.325000,95.095000,95.865000,96.635000,97.405000,98.175000,98.945000} ; is_set = true;









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
