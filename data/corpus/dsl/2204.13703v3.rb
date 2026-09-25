### Dataset description ###
# Paper: Evidence for a neutral near-threshold structure in the K_S0 recoil-mass spectra
# in e+e- -> K_S0 D_s+ D*- and e+e- -> K_S0 D_s*+ D-
# arXiv: 2204.13703v3
# Five energy points: 4.628, 4.641, 4.661, 4.682, 4.699 GeV, 3.8 fb^{-1} total

data_4628 = DatasetManager.real_data.find("705_4628")
data_4641 = DatasetManager.real_data.find("705_4641")
data_4661 = DatasetManager.real_data.find("705_4661")
data_4682 = DatasetManager.real_data.find("705_4682")
data_4699 = DatasetManager.real_data.find("705_4699")
all_data = [data_4628, data_4641, data_4661, data_4682, data_4699]

incMC_4628 = DatasetManager.inclusive_mc.find("705_4628")
incMC_4641 = DatasetManager.inclusive_mc.find("705_4641")
incMC_4661 = DatasetManager.inclusive_mc.find("705_4661")
incMC_4682 = DatasetManager.inclusive_mc.find("705_4682")
incMC_4699 = DatasetManager.inclusive_mc.find("705_4699")
all_incMC = [incMC_4628, incMC_4641, incMC_4661, incMC_4682, incMC_4699]

# Two partial reconstruction methods:
# D_s+ -tag: reconstruct bachelor K_S0 + D_s+, miss D*- -> D- pi0/gamma
# D- -tag: reconstruct bachelor K_S0 + D-, miss D_s*+ -> D_s+ gamma
# Process: e+e- -> K_S0 Z_cs(3985)^0 -> K_S0 (D_s+ D*- + D_s*+ D-)
#
# D_s+ decay modes for tag: K+K-pi+, K_S0 K+, K+K-pi+pi0, K_S0 K+pi+pi-, eta' pi+
# D- decay modes for tag: K+pi-pi-, K_S0 pi-, K_S0 pi+pi-pi-

# Decay cards for signal MC
# K_S0 D_s+ D*- channel
decay_card_dsp_dstm = <<~DECAYCARD
    Decay anti-D*-
    1.0000 anti-D0 pi- PHSP;
    Enddecay

    Decay D_s+
    1.0000 K+ K- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# K_S0 D_s*+ D- channel
decay_card_dssp_dm = <<~DECAYCARD
    Decay D_s*+
    1.0000 D_s+ gamma VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.0000 K+ K- pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_dsp_dstm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_ks_dsp_dstm"
  config.related_dataset = data_4682
  config.events = 500000
  config.decay_card = decay_card_dsp_dstm
  config.cross_section = :default
end

exMC_dssp_dm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_ks_dssp_dm"
  config.related_dataset = data_4682
  config.events = 500000
  config.decay_card = decay_card_dssp_dm
  config.cross_section = :default
end

### Event selection (BOSS) — D_s+ tag method ###
# Reconstruct bachelor K_S0 + D_s+; miss D*-
# D_s+ reconstructed via: K+K-pi+, K_S0K+, K+K-pi+pi0, K_S0K+pi+pi-, eta'pi+(eta'->pi+pi-eta, eta->gamma gamma)
alg_name_ds = "KsDsTagPartial"
alg_ds = Algorithm.new(alg_name_ds)
alg_ds.set_header(["#{alg_name_ds}Alg/#{alg_name_ds}.h"])
      .set_constant({ "ECMS" => [:double, 4.682] })

selection_ds = Selection.new
selection_ds.select_track {
              cos_theta 0.93
              Vz 10.0
              Vr 1.0
              nTot ">=4"       # K_S0 pions + D_s+ daughters
            }
            .select_photon {
              tdc_emc_start 0
              tdc_emc_end 14
              angle_to_track 10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam ">=0"
            }
            .pid(method: :chi2_sum) {
              identify :kp, :km, :pip, :pim
              chi_min_cut 0.001
            }
            .note(:ks_reconstruction, "K_S0 -> pi+ pi-: vertex fit chi2 < 100; second vertex fit (pointing to IP) chi2 < 40; decay length > 2 sigma; |M(pi+pi-) - M_K_S0| < 11 MeV/c2")
            .note(:pi0_eta_reconstruction, "pi0/eta -> gamma gamma: 1C kinematic fit constraining M(gamma gamma), chi2 < 10")
            .note(:d_ds_reconstruction, "D_s+ or D- candidates: invariant mass window within 15 MeV of known mass; best candidate chosen by mass closest to known value")
            .note(:resonance_requirements, "Table I cuts: M(K+K-)<1.05 GeV (phi), |M(K+pi-)-M(K*(892))|<70 MeV, |M(pi-pi0)-M(rho)|<150 MeV; K_S0 veto for D- -> K_S0 pi+ pi- pi-: M(pi+pi-) not in (0.48,0.52) GeV")
            .note(:recoil_mass, "Recoil mass RM(K_S0 D) = ||p_e+e- - p_K_S0_D||; resolution improved via RQ(K_S0 D) = RM(K_S0 D) + M(D) - m(D)")
            .note(:signal_region, "|RQ(K_S0 D_s+) - m(D*-)| < 20 MeV; |RQ(K_S0 D-) - m(D_s*+)| < 10 MeV")
            .note(:sideband, "M(D_s+) sideband: (1.895,1.935) + (1.995,2.035) GeV; M(D-) sideband: (1.800,1.840) + (1.900,1.940) GeV")
            .note(:partial_reconstruction, "Partial reconstruction: only bachelor K_S0 and D_s+ (or D-) are detected; the other final state particles (D*- or D_s*+ with their daughters) are not reconstructed")
            .note(:zcs_fit, "Z_cs(3985)^0 signal extracted via simultaneous unbinned maximum likelihood fit to RM(K_S0) at 5 energy points; S-wave Breit-Wigner functions used; significance 4.6 sigma")

alg_ds.with_decay_card(decay_card_dsp_dstm).apply(selection_ds)
root_files_ds = alg_ds.execute_on(all_data + all_incMC + [exMC_dsp_dstm])

### Event selection (BOSS) — D- tag method ###
# Reconstruct bachelor K_S0 + D-; miss D_s*+
# D- reconstructed via: K+pi-pi-, K_S0pi-, K_S0pi+pi-pi-
alg_name_dm = "KsDmTagPartial"
alg_dm = Algorithm.new(alg_name_dm)
alg_dm.set_header(["#{alg_name_dm}Alg/#{alg_name_dm}.h"])
      .set_constant({ "ECMS" => [:double, 4.682] })

selection_dm = Selection.new
selection_dm.select_track {
              cos_theta 0.93
              Vz 10.0
              Vr 1.0
              nTot ">=4"
            }
            .select_photon {
              tdc_emc_start 0
              tdc_emc_end 14
              angle_to_track 10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam ">=0"
            }
            .pid(method: :chi2_sum) {
              identify :kp, :km, :pip, :pim
              chi_min_cut 0.001
            }
            .note(:ks_reconstruction, "K_S0 -> pi+ pi-: vertex fit chi2 < 100; second vertex fit (pointing to IP) chi2 < 40; decay length > 2 sigma; |M(pi+pi-) - M_K_S0| < 11 MeV/c2")
            .note(:d_ds_reconstruction, "D_s+ or D- candidates: invariant mass window within 15 MeV of known mass; best candidate chosen by mass closest to known value")
            .note(:resonance_requirements, "Table I cuts: M(K+K-)<1.05 GeV (phi), |M(K+pi-)-M(K*(892))|<70 MeV, |M(pi-pi0)-M(rho)|<150 MeV/c2; K_S0 veto as above")
            .note(:recoil_mass, "Recoil mass RM(K_S0 D) = ||p_e+e- - p_K_S0_D||; RQ(K_S0 D) = RM(K_S0 D) + M(D) - m(D)")
            .note(:signal_region, "|RQ(K_S0 D-) - m(D_s*+)| < 10 MeV")
            .note(:sideband, "M(D-) sideband: (1.800,1.840) + (1.900,1.940) GeV")

alg_dm.with_decay_card(decay_card_dssp_dm).apply(selection_dm)
root_files_dm = alg_dm.execute_on(all_data + all_incMC + [exMC_dssp_dm])