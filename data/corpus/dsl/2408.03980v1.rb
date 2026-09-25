### Dataset description ###
# Ds data at multiple energies: 4.128-4.226 GeV
# Group I:   4.128 + 4.157 GeV (BOSS 705)
# Group II:  4.178 GeV (BOSS 703)
# Group III: 4.189, 4.199, 4.209, 4.219 GeV (BOSS 703)
# Group IV:  4.226 GeV (BOSS 703)
ds_data_4128 = DatasetManager.real_data.find("705_4130")     # 401.5 pb^-1
ds_data_4157 = DatasetManager.real_data.find("705_4160")     # 408.7 pb^-1
ds_data_4178 = DatasetManager.real_data.find("703_4180")     # 3189.0 pb^-1
ds_data_4189 = DatasetManager.real_data.find("703_4190")     # ~570 pb^-1
ds_data_4199 = DatasetManager.real_data.find("703_4200")     # ~526 pb^-1
ds_data_4209 = DatasetManager.real_data.find("703_4210")     # ~572 pb^-1
ds_data_4219 = DatasetManager.real_data.find("703_4220")     # ~569 pb^-1
ds_data_4226 = DatasetManager.real_data.find("703_4230")     # ~1101 pb^-1

ds_data_all = [ds_data_4128, ds_data_4157, ds_data_4178, ds_data_4189, ds_data_4199, ds_data_4209, ds_data_4219, ds_data_4226]

ds_incMC_4128 = DatasetManager.inclusive_mc.find("705_4130")
ds_incMC_4157 = DatasetManager.inclusive_mc.find("705_4160")
ds_incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
ds_incMC_4189 = DatasetManager.inclusive_mc.find("703_4190")
ds_incMC_4199 = DatasetManager.inclusive_mc.find("703_4200")
ds_incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
ds_incMC_4219 = DatasetManager.inclusive_mc.find("703_4220")
ds_incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

ds_incMC_all = [ds_incMC_4128, ds_incMC_4157, ds_incMC_4178, ds_incMC_4189, ds_incMC_4199, ds_incMC_4209, ds_incMC_4219, ds_incMC_4226]

# ------------------------------------------------------------------------------
# Decay card
# ------------------------------------------------------------------------------
decay_card_ds_signal = <<~DECAYCARD
    Decay D_s+
    1.0000  gamma   rho+                    HELAMP 1.0 0.0 0.0 0.0 0.0;
    Enddecay

    Decay rho+
    1.0000  pi+   pi0                        VSS;
    Enddecay

    Decay pi0
    1.0000  gamma   gamma                   PHSP;
    Enddecay

    End
DECAYCARD

exMC_ds_signal = DatasetManager.create_exclusive_mc_for(ds_data_all) do |config|
  config.sample_name    = "exmc_ds_gamma_rhop"
  config.events          = 500_000
  config.decay_card      = decay_card_ds_signal
  config.cross_section   = :default
end

# ==============================================================================
# TagAnalysis: D_s+ -> gamma rho(770)+  (double-tag method)
# ST D_s- reconstructed in 5 hadronic modes.
# Signal: D_s+ -> gamma rho+ -> gamma pi+ pi0
# ==============================================================================
alg_ds = TagAnalysis.new("DsGammaRhoPlus")
alg_ds.set_header(["DsGammaRhoPlusAlg/DsGammaRhoPlus.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
       .set_alias({"std::vector<double>" => "Vdouble"})

alg_ds.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKKPiPi0, :DstoKsKPiPi, :DstoKPiPi
end

alg_ds.signal_side do |s|
  s.photons 3..48              # radiative photon + 2 photons from pi0 -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(pip: 1)            # pi+ from rho+ -> pi+ pi0
  s.require_charge 1
end

alg_ds.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_ds
  .note(:recoil_mass,
        "Recoil mass M_rec against ST D_s- candidate used to select D_s*+ D_s- events; " \
        "M_rec window depends on E_cm (Table 1 in paper), applied in ROOT.")
  .note(:rho_mass_window,
        "pi+ pi0 invariant mass required in [0.62, 0.91] GeV/c^2 (rho(770)+ signal region) in ROOT.")
  .note(:eta_veto,
        "M(gamma gamma_h) > 0.68 GeV/c^2 to veto D_s+ -> pi+ eta background (ROOT).")
  .note(:extra_pi0_eta_veto,
        "M(gamma gamma_extra) veto in pi0 [0.115, 0.150] and eta [0.50, 0.57] GeV/c^2 windows (ROOT).")
  .note(:best_candidate,
        "Best candidate chosen by average mass (M_sig + M_tag)/2 closest to known D_s mass (ROOT).")
  .note(:signal_extraction,
        "2D unbinned ML fit on M_sig vs cos(theta_H) for signal extraction (ROOT).")

alg_ds.with_decay_card(decay_card_ds_signal).apply
root_files_ds = alg_ds.execute_on(ds_data_all + ds_incMC_all + [exMC_ds_signal])