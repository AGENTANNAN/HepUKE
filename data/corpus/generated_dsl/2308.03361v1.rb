# =====================================================================
# BOSS-side DSL for  e+e- -> Lambda Sigma0 + c.c.
# Measured at the 14 BESIII R-scan points spanning 2.3094 - 3.0800 GeV.
# Scope: dataset preparation + event selection (up to / instead of the
# kinematic fit).  No kinematic fit is applied anywhere - at the 13
# higher-energy points the Sigma0 is inferred from the recoil four-momentum.
# =====================================================================

### ------------------------------------------------------------------
### Dataset preparation
### ------------------------------------------------------------------
# The 14 R-scan points (BOSS 713) covering 2.3094 - 3.0800 GeV
scan_names = %w[
  713_2309 713_2386 713_2396 713_2500 713_2644 713_2646 713_2700
  713_2800 713_2900 713_2950 713_2981 713_3000 713_3020 713_3080
]

data_points  = scan_names.map { |n| DatasetManager.real_data.find(n) }        # real data
incmc_points = scan_names.map { |n| DatasetManager.inclusive_mc.find(n) }    # matching inclusive MC

# Near-threshold point (indirect method) vs. the other 13 points (single-Lambda tag)
data_2309  = data_points[0]
incmc_2309 = incmc_points[0]
data_rest  = data_points[1..]
incmc_rest = incmc_points[1..]

# ConExc decay card: models ISR up to second order and the measured Born cross
# section of e+e- -> Lambda Sigma0 (+ c.c.), with Lambda -> p pi-,
# Sigma0 -> gamma Lambda and the charge-conjugate decays.
# The DSL auto-detects the "ConExc" token, switches to the no-KKMC template and
# injects "Particle vpho <ECMS> 0.0" for every energy point, so no
# Particle vpho statement is written here.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 ConExc -1 Lambda0 Sigma0;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay Sigma0
    1.000 gamma Lambda0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.000 gamma anti-Lambda0 PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive ConExc MC for each of the 14 energy points
exmc_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "lsigma0_conexc_mc"   # auto-suffixed per energy point
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### ------------------------------------------------------------------
### Event selection - single-Lambda tag (the 13 higher-energy points)
### ------------------------------------------------------------------
alg_tag_name = "LambdaSigma0Tag"
alg_tag = Algorithm.new(alg_tag_name)
alg_tag.set_header(["#{alg_tag_name}Alg/#{alg_tag_name}.h"])
       .set_constant({ "ECMS" => [:double, 3.080] })   # placeholder - ECMS differs per scan point
       .set_alias({ "std::vector<double>" => "Vdouble" })

sel_tag = Selection.new
sel_tag.select_track {                                  # charged track selection
          cos_theta 0.93                                # |cos(theta)| < 0.93
          Vz        30.0                                # |Vz| < 30 cm
          Vr        10.0                                # Vr < 10 cm
          nChrp     ">=1"                               # at least one positive track
          nChrn     ">=1"                               # at least one negative track
        }
       .pid(method: :probability) {                     # probability PID, prob > 0.001
          prob_cut 0.001
          identify :proton, against: [:kaon, :pion]     # exactly one (anti-)proton
          identify :pim,    against: [:kaon]            # exactly one pi-
          nprp "==1"
          npim "==1"
        }
       .secondary_vertex_fit([:prp, :pim]) {            # Lambda -> p pi- (secondary vertex)
          build_virtual_particle(:Lambda).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
       .partial_rec([1])                                # reconstruct Lambda (recID 1);
                                                        # Sigma0 = P4_cms - p4_Lambda (no kinematic fit)

alg_tag
  .note(:lambda_vertex_selection, "the Lambda -> p pi- secondary-vertex fit requires decay length > 2 sigma_vtx and summed chi2 < 50; the Lambda mass window [1.11, 1.12] GeV/c^2 is imposed after by_minimizing_mass_difference")
  .note(:efficiency_curve, "the centre-of-mass energy (and hence the recoil-mass resolution) differs at every scan point; the ECMS constant below is a placeholder and must be set per dataset at run time")
  .with_decay_card(decay_card_signal).apply(sel_tag)

alg_tag.execute_on(data_rest + incmc_rest + exmc_signal[1..])

### ------------------------------------------------------------------
### Event selection - indirect method (near-threshold 2.3094 GeV)
### ------------------------------------------------------------------
alg_thr_name = "LambdaSigma0Indirect"
alg_thr = Algorithm.new(alg_thr_name)
alg_thr.set_header(["#{alg_thr_name}Alg/#{alg_thr_name}.h"])
       .set_constant({ "ECMS" => [:double, 2.3094] })   # single energy point
       .set_alias({ "std::vector<double>" => "Vdouble" })

sel_thr = Selection.new
sel_thr.select_track {                                  # charged track selection
          cos_theta 0.93                                # |cos(theta)| < 0.93
          Vz        10.0                                # |Vz| < 10 cm
          Vr        1.0                                 # Vr < 1 cm
          nChrp     ">=1"                               # at least one positive track
          nChrn     ">=1"                               # at least one negative track
        }
       .pid(method: :probability) {                     # probability PID, prob > 0.001
          prob_cut 0.001
          identify :pion, against: [:kaon, :proton]     # exactly one pi+ and one pi-
          npip "==1"
          npim "==1"
        }

alg_thr
  .note(:pid_correction_method, "the near-threshold PID uses dE/dx information only (no TOF / MUC), acceptance probability > 0.001")
  .note(:two_pion_vertex, "a two-pion (pi+ pi-) vertex is required within 2 cm of the beam line; this vertex-distance criterion has no DSL vertex-block equivalent")
  .note(:soft_pion_momentum, "the pi- momentum is required to lie in 0.08-0.12 GeV/c")
  .note(:antiproton_annihilation, "at least two additional secondary tracks must originate from a common vertex 1-5 cm from the beam line (antiproton-annihilation signature)")
  .note(:no_kinematic_fit, "indirect method: no kinematic fit is applied; the Lambda Sigma0 signal is extracted from a fit to the pi+ momentum spectrum in the ROOT stage")
  .with_decay_card(decay_card_signal).apply(sel_thr)

alg_thr.execute_on([data_2309, incmc_2309, exmc_signal[0]])