# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Data at E_cm = 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV (7.33 fb^-1 total).
ds_energies = ["703_4128", "703_4157", "703_4178", "703_4189", "703_4199", "703_4209", "703_4219", "703_4226"]
ds_datasets = ds_energies.map { |e| DatasetManager.real_data.find(e) }
ds_incMCs   = ds_energies.map { |e| DatasetManager.inclusive_mc.find(e) }

# Signal decay cards for the two K*(892)+ sub-modes
decay_card_kpi0 = <<~DECAYCARD
    Decay D_s*+
    1.0000 gamma D_s+                       VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.0000 gamma K*+                        SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay K*+
    1.0000 K+ pi0                           VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                      PHSP;
    Enddecay

    Decay D_s-
    1.0000 anti-K0 K-                       PHSP;
    Enddecay

    End
DECAYCARD

decay_card_kspi = <<~DECAYCARD
    Decay D_s*+
    1.0000 gamma D_s+                       VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.0000 gamma K*+                        SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay K*+
    1.0000 K_S0 pi+                         VSS;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                          PHSP;
    Enddecay

    Decay D_s-
    1.0000 anti-K0 K-                       PHSP;
    Enddecay

    End
DECAYCARD

exMC_kpi0 = DatasetManager.create_exclusive_mc_for(ds_datasets) do |config|
  config.sample_name    = "Ds_gamma_Kstar_KPi0"
  config.events         = 150000
  config.decay_card     = decay_card_kpi0
  config.cross_section  = :default
end
exMC_kpi0.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

exMC_kspi = DatasetManager.create_exclusive_mc_for(ds_datasets) do |config|
  config.sample_name    = "Ds_gamma_Kstar_KsPi"
  config.events         = 150000
  config.decay_card     = decay_card_kspi
  config.cross_section  = :default
end
exMC_kspi.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Tag Analysis 1: K*+ -> K+ pi0 ###
alg_kpi0 = TagAnalysis.new("DsTagRadKKpi0")
alg_kpi0.set_header(["DsTagRadKKpi0Alg/DsTagRadKKpi0.h"])
        .set_constant({"ECMS" => [:double, 4.178]})
        .set_alias({"std::vector<double>" => "Vdouble"})

# Single-tag Ds- with three hadronic modes
alg_kpi0.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPiPi
  t.charm -1
end

# Signal side: gamma (radiative), K+ from K*+ -> K+ pi0, and pi0 -> gamma gamma
alg_kpi0.signal_side do |s|
  s.photons 3                       # >=3 photons: the radiative gamma + 2 from pi0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(kp: 1)                  # signal side: one K+ from K*+ -> K+ pi0
end

alg_kpi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_kpi0.note(:radiative_photon,
              "The radiative photon is taken to be the most energetic photon in the event; require E_gamma > 0.55 GeV")
        .note(:kstar_window,
              "Require M(K+ pi0) in (0.83, 0.94) GeV/c^2 (K*(892)+ window)")
        .note(:extra_photon_vetoes,
              "Veto M(gamma gamma_extra) in the pi0 [0.115, 0.150] and eta [0.50, 0.57] GeV/c^2 windows; require M(gamma gamma_h) > 0.62 GeV/c^2 to reject Ds+ -> K+ eta background")
        .note(:tag_recoil_window,
              "Retain events with M_rec against ST Ds- within the per-energy window given in Table 1 (2.040-2.220 GeV/c^2 depending on E_cm)")

alg_kpi0.apply
alg_kpi0.execute_on(ds_datasets + ds_incMCs + exMC_kpi0)

### Tag Analysis 2: K*+ -> K_S0 pi+ ###
alg_kspi = TagAnalysis.new("DsTagRadKsPi")
alg_kspi.set_header(["DsTagRadKsPiAlg/DsTagRadKsPi.h"])
        .set_constant({"ECMS" => [:double, 4.178]})
        .set_alias({"std::vector<double>" => "Vdouble"})

alg_kspi.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPiPi
  t.charm -1
end

# Signal side: gamma (radiative) + pi+ + K_S0(-> pi+ pi-)
alg_kspi.signal_side do |s|
  s.photons 1                       # radiative gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(pip: 2, pim: 1)         # pi+ from K*+ + 2 pions from K_S0
end

alg_kspi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

alg_kspi.note(:radiative_photon,
              "Radiative gamma is the most energetic photon; E_gamma > 0.55 GeV")
        .note(:kstar_window,
              "Require M(K_S0 pi+) in (0.83, 0.94) GeV/c^2 for the K*+")
        .note(:ks_reconstruction,
              "K_S0 reconstructed from two oppositely charged tracks (assigned as pi+ pi-, no PID) with |Vz| < 20 cm, vertex-fit chi^2 < 100, M(pi+pi-) in [0.487, 0.511] GeV/c^2, and decay length > 2 sigma from IP")
        .note(:extra_photon_vetoes,
              "Veto M(gamma gamma_extra) in the pi0 [0.115, 0.150] and eta [0.50, 0.57] GeV/c^2 windows")
        .note(:tag_recoil_window,
              "Retain events with M_rec against ST Ds- within the per-energy window given in Table 1")

alg_kspi.apply
alg_kspi.execute_on(ds_datasets + ds_incMCs + exMC_kspi)
