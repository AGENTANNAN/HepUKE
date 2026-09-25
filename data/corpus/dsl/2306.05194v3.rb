# Paper: 2306.05194v3
# Title: Precision Measurements of Ds+ -> eta e+ nu_e and Ds+ -> eta' e+ nu_e
# Energy: 4.128-4.226 GeV (8 energy points)
# Double-tag (DT) method:
#   ST: Ds- -> 14 hadronic tag modes
#   DT: Ds+ -> eta(') e+ nu_e (semileptonic signal)
#   eta -> gamma gamma or pi0 pi+ pi-
#   eta' -> eta pi+ pi- or gamma rho0
# Semileptonic: missing neutrino

### Dataset preparation ###
data_703_4178 = DatasetManager.load_real_data.find("703_4178")

all_data = [data_703_4178]
all_incMC = DatasetManager.load_inclusive_mc

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Ds_to_eta_etap_e_nu"
  config.events         = 500000
  config.decay_card     = nil  # TagAnalysis
end

### TagAnalysis: Double-tag with Ds- ST hadronic tags ###
alg = TagAnalysis.new("Ds_to_eta_etap_e_nu")
alg.set_header(["DsEtaEtapENuAlg/DsEtaEtapENu.h"])

# ST: Ds- -> 14 hadronic tag modes
# Representative subset of well-known Ds tag modes
alg.tag_side(:Dsm) {
  modes(:DstoKKPi, :DstoKPiPi, :DstoPiPiPi,
        :DstoKKPiPi0, :DstoKSPiPiPi)
  charm -1
}

# Signal side: Ds+ -> eta(') e+ nu_e
# eta -> gamma gamma; eta' -> eta pi+ pi- or gamma rho0
# Missing neutrino; e+ identified via PID
# Transition gamma/pi0 from Ds*+ -> gamma/pi0 Ds+
alg.signal_side {
  charged(ep: 1)
  photons 3     # eta -> gamma gamma plus transition gamma
  missing :nu_e
}

# Fit: 3C kinematic fit: energy-momentum conservation + constrain
# both Ds masses + constrain Ds* mass
# Then M_miss^2 fit for signal yield extraction
alg.fit {
  constrain_four_momentum
  chi2_cut 200
}

# ST: Ds- hadronic tags, M_BC + M_tag fit
# 14 tag modes with M_BC and M_tag signal regions
alg.note(:st_selection,
  "ST: 14 Ds- hadronic tag modes. Track: |cos theta|<0.93, Vz<10 cm, Vxy<1 cm. PID: L_K>L_pi for kaons, L_pi>L_K for pions. K_S0: pi+pi- vertex fit, |M-M(K_S0)|<12 MeV/c2. pi0/eta: gamma gamma, mass windows. M_BC cuts per energy point. M_tag fits for ST yields. Applied in ROOT.")

# DT: semileptonic signal, transition gamma/pi0, kinematic fit,
# M_miss^2 simultaneous fit for eta and eta' channels
alg.note(:dt_selection,
  "DT: Ds+ -> eta(') e+ nu_e. e+ PID: CL_e>0.001, CL_e/(CL_e+CL_pi+CL_K)>0.8. Transition gamma/pi0 with smallest |DeltaE|. 3C kinematic fit: energy-momentum conservation + Ds/Ds* mass constraints. chi2<200 for eta'_gamma_rho0 mode. Extra energy and track vetoes. M(eta(')e+)<1.9 GeV/c2. M_miss^2 simultaneous fit. BF(eta e nu) = (2.255+/-0.039+/-0.051)%. BF(eta' e nu) = (0.810+/-0.038+/-0.024)%. Applied in ROOT.")

# Form factor extraction and eta-eta' mixing angle
alg.note(:form_factors,
  "Hadronic transition form factors extracted via fits to partial decay rates in q^2 bins. Modified pole model and 2-Par series expansion. Products f_+^eta(0)|V_cs| and f_+^eta'(0)|V_cs| determined. eta-eta' mixing angle phi_P = (40.0+/-2.0+/-0.6) degrees. Applied in ROOT.")

# 7.33 fb^-1 total at 4.128-4.226 GeV
alg.note(:energy_points,
  "8 c.m. energies: 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV. Total 7.33 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on([all_datasets])