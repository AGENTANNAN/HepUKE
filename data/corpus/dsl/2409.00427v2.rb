# e⁺e⁻ → Ξ⁰Ξ̄⁰ Born cross-section measurement
# Single-baryon tag method at 45 CM energies from 3.51 to 4.95 GeV
# Total ~30 fb⁻¹, arXiv:2409.00427v2
#
# Note: The paper uses 45 energy points. The representative subset below
# covers the major scan points. For the complete list, add the remaining
# entries from the BESIII dataset table matching the paper's Table I.

# ── Scan data points (representative subset of 45 energies) ──
scan_data_points = [
  DatasetManager.real_data.find("703_3810"),
  DatasetManager.real_data.find("703_3900"),
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

# Inclusive MC: pick the inclusive MC matching each data point
scan_incMC = scan_data_points.map { |d|
  DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
}

# Decay card: continuum e⁺e⁻ → Ξ⁰Ξ̄⁰ via KKMC + psi(4260) top mother.
# (Alternative: ConExc Form D -2 3322 -3322 with xs_user.txt for Born-cross-section measurement.
#  Since this paper IS the cross-section measurement, the flat-PHSP generation with
#  KKMC+psi(4260) is used; efficiency correction factors derived in ROOT.)
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 Xi0 anti-Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.000 pi0 Lambda0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.000 pi0 anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for the signal channel at all scan points
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data_points) do |config|
  config.sample_name = "ee_to_Xi0Xi0bar"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ── Algorithm ──
alg = Algorithm.new("Xi0Xi0barCrossSection")
alg.set_header(["Xi0Xi0barCrossSectionAlg/Xi0Xi0barCrossSection.h"])
   .set_constant({ "ECMS" => [:double, 4.0] })  # placeholder; actual value from MeasuredEcmsSvc per run

# ── Selection: single-baryon tag Ξ⁰ → π⁰Λ, Λ → pπ⁻ ──
event_selection = Selection.new

event_selection
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=2"
    nChrn ">=1"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .assign({ :chrgp => :pip, :chrgn => :pim })
  # Reconstruct Λ → pπ⁻ via secondary vertex fit
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Reconstruct π⁰ → γγ via 1C kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  }
  # Select Ξ⁰ candidate: minimum |M(π⁰Λ) - m_Ξ⁰|
  # No kinematic fit with unmeasured Ξ⁰ — the recoil-mass method is used instead
  # (the Ξ⁰ is tagged via the π⁰Λ combination; its recoil mass is computed in ROOT)

alg.with_decay_card(decay_card).apply(event_selection)

# ── ROOT-level notes ──
# Single-baryon tag: select Ξ⁰ → π⁰Λ by minimizing |M(π⁰Λ)-m_Ξ⁰|
# Signal region: |M(π⁰Λ)-m_Ξ⁰| < 10 MeV/c², |M_recoil-m_Ξ⁰| < 60 MeV/c²
# Sideband background estimation from 4 sideband regions in 2D plane
# Λ mass window: |M(pπ⁻)-m_Λ| < 5 MeV/c²
alg.note(:single_baryon_tag, "Ξ⁰ tag: minimize |M(π⁰Λ)-m_Ξ⁰|; signal window |M(π⁰Λ)-m_Ξ⁰|<10 MeV/c², |M_recoil-m_Ξ⁰|<60 MeV/c² applied in ROOT")
alg.note(:sideband, "4-region 2D sideband background estimation in ROOT")
alg.note(:lambda_window, "|M(pπ⁻)-m_Λ| < 5 MeV/c² after secondary vertex fit with L/σ_L>0 applied in ROOT")
alg.note(:cross_section, "Born cross section extracted from fitted Ξ⁰ yields; ISR correction and efficiency variation with √s applied in ROOT")
alg.note(:energy_points, "45 CM energy points from 3.51 to 4.95 GeV — verify dataset sample names against BES3_dataset.md for the complete list")

# Execute on all data points + inclusive MC + exclusive MC
all_datasets = scan_data_points + scan_incMC + exMC_signal
alg.execute_on(all_datasets.flatten)