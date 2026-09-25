# DSL for paper 2406.01332v3: Measurements of the branching fractions of semileptonic
# Ds+ decays via e+e- → Ds*+Ds*- using a double-tag method
# 14 ST tag modes for Ds*-; 6 semileptonic signal modes: ηe+ν, η'e+ν, φe+ν, f0e+ν, K0e+ν, K*0e+ν
# 10.64 fb-1 at √s = 4.237-4.699 GeV

# Data: BOSS 703/705 scan points
scan_data = DatasetManager.real_data.where(cms_energy: {value: 4237..4700})
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: {value: 4237..4700})

# Decay card for signal MC: e+e- → Ds*+Ds*- via KKMC
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s*- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 K+ K- pi+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC for each energy scan point
scan_exMC = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_DsDsStarStar"
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

# TagAnalysis: double-tag Ds*- → γ(π0)Ds-, Ds- in hadronic modes
# Signal side: Ds+ → hadron e+ νe (semileptonic)
alg = TagAnalysis.new("DsDsSemileptonic")
alg.set_header(["DsDsSemileptonicAlg/DsDsSemileptonic.h"])
   .set_constant({ "ECMS" => [:double, 4.400] })
   .note(:double_tag_method, "Double-tag (DT) method: fully reconstruct Ds*- in hadronic tag mode (ST), then reconstruct semileptonic Ds+ decay on signal side; branching fraction B_sig = N_DT / (N_ST × ε_sig × B_sub)")
   .note(:tag_modes, "14 ST Ds*- tag modes: Ds*- → γ(π0)Ds-, with Ds- → K+K-π-, K+K-π-π0, π+π-π-, KS0K-, KS0K-π0, KS0KS0π-, KS0K+π-π-, KS0K-π+π-, K-π+π-, ηγγπ-, ηπ+π-π0π-, η'π+π-ηπ-, η'γρ0π-, ηγγρ-; ST Ds*+ reconstructed via Ds*+ → γ(π0)Ds+")
   .note(:st_selection, "ST candidates selected via ΔE and M_BC requirements; M_BC signal region varies with energy (2.104-2.123 GeV/c²); best candidate per tag mode chosen by minimum |ΔE|; ST yields from fits to M_BC spectra")
   .note(:signal_side, "Signal side uses surviving tracks/showers not used in ST; Ds+ → ηe+νe, η'e+νe, φe+νe, f0(980)e+νe, K0e+νe, K*0e+νe; e+ PID: CL_e > 0.001 and CL_e/(CL_e+CL_π+CL_K) > 0.8; bremsstrahlung recovery by adding EMC showers within 10°")
   .note(:signal_side_resonances, "Intermediate resonances: φ → K+K- mass [1.004, 1.034], f0(980) → π+π- mass [0.880, 1.080], K*0 → K+π- mass [0.882, 0.992] GeV/c²")
   .note(:missing_mass, "Signal yield extracted from fit to M²_miss distribution; M²_miss peaks at zero for signal; Ds*+ momentum constrained using ST Ds*- direction; invariant mass of hadron+lepton < 1.90 (1.75) GeV/c² for Cabibbo-favored (suppressed) decays")
   .note(:extra_energy_veto, "Maximum unused shower energy E_max_extra_γ < 0.3 GeV; no extra charged tracks allowed on DT side")
   .note(:simultaneous_fit, "For Ds+ → ηe+νe and η'e+νe, simultaneous fits to two η/η' decay modes sharing common branching fraction; ROOT level")
   .note(:hadronic_form_factors, "Hadronic form factors determined via two-parameter series expansion fit to partial decay rates in q² intervals; ROOT level")
   .with_decay_card(decay_card)

# Tag side: single tag Ds*-
alg.tag_side(:Ds_star) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi,
          :DstoKsK, :DstoKsKPi0, :DstoKsKsPi,
          :DstoKsKPiPi, :DstoKsKPiPi2, :DstoKPiPi,
          :DstoEtaPi, :DstoEtaPiPi0Pi, :DstoEtapEtaPi,
          :DstoEtapRhoPi, :DstoEtaRho
  t.charm -1
end

# Signal side: semileptonic Ds+ → hadron e+ νe
alg.signal_side do |s|
  s.charged :ep
  s.missing :nu_e
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply
alg.execute_on([scan_data, scan_incMC, scan_exMC])