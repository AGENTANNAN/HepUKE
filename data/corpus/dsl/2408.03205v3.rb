# DSL for BESIII paper 2408.03205v3: Sigma+ transverse polarization measurement
# e+e- -> gamma*/Psi -> Sigma+ Sigma- -> p anti-p pi0 pi0
# Seven energy points: 3.682, 3.683, 3.684, 3.685, 3.687, 3.691, 3.710 GeV
# Total integrated luminosity: 652.1 pb-1

### Dataset preparation ###
# Main datasets around psi(3686) resonance
ds_3682 = DatasetManager.real_data.find("709_3682")
ds_3683 = DatasetManager.real_data.find("709_3683")
ds_3684 = DatasetManager.real_data.find("709_3684")
ds_3685 = DatasetManager.real_data.find("709_3685")
ds_3687 = DatasetManager.real_data.find("709_3687")
ds_3691 = DatasetManager.real_data.find("709_3691")
ds_3710 = DatasetManager.real_data.find("709_3710")

scan_datasets = [ds_3682, ds_3683, ds_3684, ds_3685, ds_3687, ds_3691, ds_3710]

decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Sigma+ anti-Sigma-           PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0                       PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0                  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the energy scan
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name = "SigmaPlus_SigmaMinus_p_antip_pi0_pi0"
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("SigmaPlusTransPolar")
alg.set_header(["SigmaPlusTransPolarAlg/SigmaPlusTransPolar.h"])

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  nChrp ">=1"
                  nChrn ">=1"
                  nNet "==0"
                }
               .select_photon {
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  tdc_emc_start 0
                  tdc_emc_end 14
                  nGam ">=4"
                }
               # No PID required: (anti-)protons from Sigma+ decay identified by momentum > 0.5 GeV/c
               .assign({chrgp: :prp, chrgn: :prm})
               # Reconstruct two pi0 from photon pairs
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 "==2"
                }
               # 6C kinematic fit: 4-momentum + 2 pi0 mass constraints
               .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 100
                }

# Set ECMS dynamically via a constant (single algorithm handles all energies)
alg.set_constant({"ECMS" => [:double, 3.686]})

alg.note(:no_vertex_requirement, "no Vz/Vr cut applied to avoid efficiency loss")
   .note(:proton_momentum_cut, "proton/anti-proton momentum > 0.5 GeV/c instead of PID")
   .note(:best_combination, "combination with smallest chi2_6C retained; paired by min sqrt((M(ppi0)-mSigma)^2 + (M(anti-p pi0)-mAntiSigma)^2)")
   .note(:signal_region, "M(ppi0) and M(anti-p pi0) within [mSigma-4sigma, mSigma+3sigma]")
   .note(:jpsi_veto, "|M(pi0 pi0)_recoil - m(J/psi)| > 15 MeV/c^2 to suppress e+e- -> pi0 pi0 J/psi")
   .with_decay_card(decay_card)
   .apply(event_selection)

root_files = alg.execute_on(scan_datasets.map { |ds| [ds] }.flatten)
# Note: inclusive MC and exclusive MC datasets also to be added:
# scan_incMCs = scan_datasets.map { |ds| DatasetManager.inclusive_mc.find(ds.sample_name) }
# root_files = alg.execute_on(scan_datasets + scan_incMCs + exMCs_signal)