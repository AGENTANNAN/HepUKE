### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Matching inclusive MC sample

# Decay card for the signal process ψ(2S) → γ ηc , ηc → 2(π+π-) η , η → γγ (phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 pi+ pi- pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_gamma_etac_2pipi_eta_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "PsipGamma4PiEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                # Charged-track selection
                  cos_theta 0.93              # |cosθ| < 0.93
                  Vz        10.0              # |Vz| < 10 cm
                  Vr        1.0               # Vr < 1 cm in the transverse plane
                  nChrp     "==2"             # exactly two π+ candidates
                  nChrn     "==2"             # exactly two π- candidates
                  nNet      "==0"             # net charge zero
                }
               .select_photon {               # Photon selection
                  tdc_emc_start     0         # EMC TDC start time
                  tdc_emc_end       14        # EMC TDC end time
                  angle_to_track    10.0      # ≥ 10° separation from any charged track
                  energyThreshold_b 0.025     # 25 MeV in the barrel
                  energyThreshold_e 0.050     # 50 MeV in the endcap
                  nGam              ">=3"     # at least three photons (1 radiative γ + 2 from η)
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit to reconstruct η → γγ
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 20                 # χ² < 20 for the mass-constrained fit
                  neta     ">=1"              # require at least one η candidate
                }
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks assumed π+/π-
               .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :eta]) {   # 4C fit of γ π+π+π-π- η
                  nominal                     # nominal fit (corrected four-momenta come from here)
                  constrain_four_momentum     # 4-momentum conservation against the CMS energy
                  chi2_cut 200                # loose χ² cut; optimal 5C χ² < 15 applied later in ROOT
                }
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {  # 4γ hypothesis
                  constrain_four_momentum     # no chi2_cut / no nominal: stores competing χ² for ROOT veto
                }
               .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {                 # 2γ hypothesis
                  constrain_four_momentum     # stores competing χ² for ROOT veto
                }

my_algorithm
  .note(:background_veto_psip_pjpsi,
        "ψ(3686)→P J/ψ vetoed: events whose recoil mass against the prompt hadron system P " \
        "lies inside the J/ψ mass window are rejected (applied as a post-kinematic-fit cut in ROOT).")
  .note(:background_veto_psip_pi0h,
        "ψ(3686)→π0 H vetoed: events where the γγ invariant mass falls in the π0 region are rejected.")
  .note(:background_veto_gamma_chicj,
        "For the γηc(2S)/χcJ modes, events consistent with ψ(3686)→γχcJ, χcJ→γJ/ψ are vetoed.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])