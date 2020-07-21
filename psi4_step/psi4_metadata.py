methods = {
    'effective fragment potential (EFP)': {
        'method': 'efp',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    'Kohn-Sham (KS) density functional theory (DFT)': {
        'method': 'dft',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'Hartree-Fock (HF) self consistent field (SCF)': {
        'method': 'hf',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'HF with dispersion, BSSE, and basis set corrections': {
        'method': 'hf3c',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'PBEh with dispersion, BSSE, and basis set corrections': {
        'method': 'pbeh3c',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'density cumulant functional theory': {
        'method': 'dcft',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    '2nd-order Møller–Plesset perturbation theory (MP2)': {
        'method': 'mp2',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    '3rd-order Møller–Plesset perturbation theory (MP3)': {
        'method': 'mp3',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'MP3 with frozen natural orbitals': {
        'method': 'fno-mp3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'average of MP2 and MP3': {
        'method': 'mp2.5',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    '4th-order MP perturbation theory (MP4) less triples': {
        'method': 'mp4(sdq)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'MP4 (less triples) with frozen natural orbitals': {
        'method': 'fno-mp4(sdq)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'full MP4': {
        'method': 'mp4',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'full MP4 with frozen natural orbitals': {
        'method': 'fno-mp4',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'nth-order Møller–Plesset (MP) perturbation theory': {
        'method': 'mpn',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'nth-order z-averaged perturbation theory (ZAPT)': {
        'method': 'zaptn',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'orbital-optimized second-order MP perturbation theory': {
        'method': 'omp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    'spin-component scaled OMP2': {
        'method': 'scs-omp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'a special version of SCS-OMP2 for nucleobase interactions': {
        'method': 'scs(n)-omp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'a special version of SCS-OMP2 (from ethene dimers)': {
        'method': 'scs-omp2-vdw',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'spin-opposite scaled OMP2': {
        'method': 'sos-omp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'A special version of SOS-OMP2 for pi systems': {
        'method': 'sos-pi-omp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'orbital-optimized third-order MP perturbation theory': {
        'method': 'omp3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    'spin-component scaled OMP3': {
        'method': 'scs-omp3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'a special version of SCS-OMP3 for nucleobase interactions': {
        'method': 'scs(n)-omp3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'a special version of SCS-OMP3 (from ethene dimers)': {
        'method': 'scs-omp3-vdw',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'spin-opposite scaled OMP3': {
        'method': 'sos-omp3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'A special version of SOS-OMP3 for pi systems': {
        'method': 'sos-pi-omp3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
        },
    'orbital-optimized MP2.5': {
        'method': 'omp2.5',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    'coupled electron pair approximation variant 0': {
        'method': 'cepa(0)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CEPA(0) with frozen natural orbitals': {
        'method': 'fno-cepa(0)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'coupled electron pair approximation variant 1': {
        'method': 'cepa(1)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CEPA(1) with frozen natural orbitals': {
        'method': 'fno-cepa(1)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'coupled electron pair approximation variant 3': {
        'method': 'cepa(3)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CEPA(3) with frozen natural orbitals': {
        'method': 'fno-cepa(3)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'averaged coupled-pair functional': {
        'method': 'acpf',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'ACPF with frozen natural orbitals': {
        'method': 'fno-acpf',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'averaged quadratic coupled cluster': {
        'method': 'aqcc',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'AQCC with frozen natural orbitals': {
        'method': 'fno-aqcc',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'quadratic CI singles doubles (QCISD)': {
        'method': 'qcisd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'QCISD with frozen natural orbitals': {
        'method': 'fno-qcisd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'Linear CCD': {
        'method': 'lccd',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'LCCD with frozen natural orbitals': {
        'method': 'fno-lccd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'orbital optimized LCCD': {
        'method': 'olccd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    'approximate coupled cluster singles and doubles (CC2)': {
        'method': 'cc2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'coupled cluster doubles (CCD)': {
        'method': 'ccd',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'coupled cluster singles and doubles (CCSD)': {
        'method': 'ccsd',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'Brueckner coupled cluster doubles (BCCD)': {
        'method': 'bccd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CCSD with frozen natural orbitals': {
        'method': 'fno-ccsd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'QCISD with perturbative triples': {
        'method': 'qcisd(t)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'QCISD(T) with frozen natural orbitals': {
        'method': 'fno-qcisd(t)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CCSD with perturbative triples (CCSD(T))': {
        'method': 'ccsd(t)',
        'calculation': ['energy', 'gradients'],
        'level': 'normal',
        'gradients': 'analytic'
    },
    'CCSD with asymmetric perturbative triples (CCSD(AT))': {
        'method': 'ccsd(at)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'BCCD with perturbative triples': {
        'method': 'bccd(t)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CCSD(T) with frozen natural orbitals': {
        'method': 'fno-ccsd(t)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'approximate CC singles, doubles, and triples (CC3)': {
        'method': 'cc3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'expert full control over ccenergy module': {
        'method': 'ccenergy',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'expert full control over dfocc module': {
        'method': 'dfocc',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'configuration interaction (CI) singles and doubles (CISD)': {
        'method': 'cisd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CISD with frozen natural orbitals': {
        'method': 'fno-cisd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CI singles, doubles, and triples (CISDT)': {
        'method': 'cisdt',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'CI singles, doubles, triples, and quadruples (CISDTQ)': {
        'method': 'cisdtq',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'nth-order CI': {
        'method': 'cin',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'full configuration interaction (FCI)': {
        'method': 'fci',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'expert full control over detci module': {
        'method': 'detci',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'complete active space self consistent field (CASSCF)': {
        'method': 'casscf',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'restricted active space self consistent field (RASSCF)': {
        'method': 'rasscf',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'multiconfigurational self consistent field (SCF)': {
        'method': 'mcscf',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'Mukherjee multireference coupled cluster (Mk-MRCC)': {
        'method': 'psimrcc',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'density matrix renormalization group SCF': {
        'method': 'dmrg-scf',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'density matrix renormalization group CASPT2': {
        'method': 'dmrg-caspt2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'density matrix renormalization group CI': {
        'method': 'dmrg-ci',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
        },
    '0th-order symmetry adapted perturbation theory (SAPT)': {
        'method': 'sapt0',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    '0th-order SAPT with special exchange scaling': {
        'method': 'ssapt0',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    '0th-order functional and/or intramolecular SAPT': {
        'method': 'fisapt0',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    '2nd-order SAPT, traditional definition': {
        'method': 'sapt2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT including all 2nd-order terms': {
        'method': 'sapt2+',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT including perturbative triples': {
        'method': 'sapt2+(3)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT including all 3rd-order terms': {
        'method': 'sapt2+3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+ with CC-based dispersion': {
        'method': 'sapt2+(ccd)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+(3) with CC-based dispersion': {
        'method': 'sapt2+(3)(ccd)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+3 with CC-based dispersion': {
        'method': 'sapt2+3(ccd)',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT including all 2nd-order terms and MP2 correction': {
        'method': 'sapt2+dmp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT including perturbative triples and MP2 correction': {
        'method': 'sapt2+(3)dmp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT including all 3rd-order terms and MP2 correction': {
        'method': 'sapt2+3dmp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+ with CC-based dispersion and MP2 correction': {
        'method': 'sapt2+(ccd)dmp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+(3) with CC-based dispersion and MP2 correction': {
        'method': 'sapt2+(3)(ccd)dmp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+3 with CC-based dispersion and MP2 correction': {
        'method': 'sapt2+3(ccd)dmp2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    '0th-order SAPT plus charge transfer (CT) calculation': {
        'method': 'sapt0-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2 plus CT': {
        'method': 'sapt2-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+ plus CT': {
        'method': 'sapt2+-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+(3) plus CT': {
        'method': 'sapt2+(3)-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+3 plus CT': {
        'method': 'sapt2+3-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+(CCD) plus CT': {
        'method': 'sapt2+(ccd)-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+(3)(CCD) plus CT': {
        'method': 'sapt2+(3)(ccd)-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'SAPT2+3(CCD) plus CT': {
        'method': 'sapt2+3(ccd)-ct',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    '2nd-order algebraic diagrammatic construction (ADC)': {
        'method': 'adc',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'EOM-CC2': {
        'method': 'eom-cc2',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    },
    'equation of motion (EOM) CCSD': {
        'method': 'eom-ccsd',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'analytic'
    },
    'EOM-CC3': {
        'method': 'eom-cc3',
        'calculation': ['energy', 'gradients'],
        'level': 'expert',
        'gradients': 'finite-difference'
    }
}  # yapf: disable

dft_functionals = {
    'B1LYP Hyb-GGA Exchange-Correlation Functional': {
        'method': 'b1lyp',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'B1PW91 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b1pw91',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B1WC Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b1wc',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B2GPPLYP Double Hybrid Exchange-Correlation Functional': {
        'name': 'b2gpplyp',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'B2PLYP Double Hybrid Exchange-Correlation Functional': {
        'name': 'b2plyp',
        'dispersion': ['none', 'd3bj', 'd3mbj', 'nl'],
        'level': 'expert',
    },
    'B3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b3lyp',
        'dispersion': ['none', 'd3bj', 'd3mbj', 'nl'],
        'level': 'normal',
    },
    'B3LYP5 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b3lyp5',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B3LYPs Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b3lyps',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B3P86 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b3p86',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'B3PW91 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b3pw91',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'B5050LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b5050lyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B86B95 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b86b95',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B86BPBE GGA Exchange-Correlation Functional': {
        'name': 'b86bpbe',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B88B95 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b88b95',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'B97-0 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b97-0',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B97-1 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b97-1',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'B97-1p Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b97-1p',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B97-2 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b97-2',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'B97-3 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b97-3',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'b97-d': {
        'name': 'b97-d',
        'dispersion': ['none', 'd3bj', 'd3mbj'],
        'level': 'expert',
    },
    'B97-GGA1 GGA Exchange-Correlation Functional': {
        'name': 'b97-gga1',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B97-K Hyb-GGA Exchange-Correlation Functional': {
        'name': 'b97-k',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'b97m-d3bj': {
        'name': 'b97m-d3bj',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'B97M-V GGA Exchange-Correlation Functional': {
        'name': 'b97m-v',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'BB1K Hyb-GGA Exchange-Correlation Functional': {
        'name': 'bb1k',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'BHandH Hyb-GGA Exchange-Correlation Functional': {
        'name': 'bhandh',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'BHandHLYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'bhandhlyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'BLYP GGA Exchange-Correlation Functional': {
        'name': 'blyp',
        'dispersion': ['none', 'd3bj', 'd3mbj', 'nl'],
        'level': 'expert',
    },
    'BOP GGA Exchange-Correlation Functional': {
        'name': 'bop',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'BP86 GGA Exchange-Correlation Functional': {
        'name': 'bp86',
        'dispersion': ['none', 'd3bj', 'd3mbj', 'nl'],
        'level': 'expert',
    },
    'CAM-B3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'cam-b3lyp',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'CAP0 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'cap0',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'core-dsd-blyp': {
        'name': 'core-dsd-blyp',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'Dispersionless Hybrid Meta-GGA XC Functional': {
        'name': 'dldf',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Dispersionless Hybrid Meta-GGA XC Functional + d09': {
        'name': 'dldf+d09',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Dispersionless Hybrid Meta-GGA XC Functional + d10': {
        'name': 'dldf+d10',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'DSD-BLYP SCS Double Hybrid XC Functional': {
        'name': 'dsd-blyp',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'DSD-PBEB95 SCS Double Hybrid Meta-GGA XC Functional': {
        'name': 'dsd-pbeb95',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'DSD-PBEP86 SCS Double Hybrid XC Functional': {
        'name': 'dsd-pbep86',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'DSD-PBEPBE SCS Double Hybrid XC Functional': {
        'name': 'dsd-pbepbe',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'EDF1 GGA Exchange-Correlation Functional': {
        'name': 'edf1',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'EDF2 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'edf2',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'FT97 GGA Exchange-Correlation Functional': {
        'name': 'ft97',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'GAM GGA Minnesota Exchange-Correlation Functional': {
        'name': 'gam',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HCTH120 GGA Exchange-Correlation Functional': {
        'name': 'hcth120',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'HCTH147 GGA Exchange-Correlation Functional': {
        'name': 'hcth147',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HCTH407 GGA Exchange-Correlation Functional': {
        'name': 'hcth407',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'HCTH407P GGA Exchange-Correlation Functional': {
        'name': 'hcth407p',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HCTH93 GGA Exchange-Correlation Functional': {
        'name': 'hcth93',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HCTHP14 GGA Exchange-Correlation Functional': {
        'name': 'hcthp14',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HCTHP76 GGA Exchange-Correlation Functional': {
        'name': 'hcthp76',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'hf': {
        'name': 'hf',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'hf+d': {
        'name': 'hf+d ',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Hartree Fock based 3C composite method with minimal basis set, gCP and D3(BJ)': {  # noqa: E501
        'name': 'hf3c',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HJS-B88 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hjs-b88',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HJS-B97X Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hjs-b97x',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HJS-PBE Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hjs-pbe',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HJS-PBE-SOL Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hjs-pbe-sol',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HPBEINT Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hpbeint',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'HSE03 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hse03',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'HSE06 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'hse06',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'KSDT Exchange-Correlation Functional': {
        'name': 'ksdt',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'KT2 GGA Exchange-Correlation Functional': {
        'name': 'kt2',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'LC-VV10 GGA Exchange-Correlation Functional': {
        'name': 'lc-vv10',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'LRC-WPBE GGA Exchange-Correlation Functional': {
        'name': 'lrc-wpbe',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'LRC-wPBEh Hyb-GGA Exchange-Correlation Functional': {
        'name': 'lrc-wpbeh',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'M05 Meta-GGA XC Functional': {
        'name': 'm05',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Heavily Parameterized Hybrid M05-2X Meta-GGA XC Functional': {
        'name': 'm05-2x',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'M06 Meta-GGA XC Functional': {
        'name': 'm06',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Hybrid M06-2X Meta-GGA XC Functional': {
        'name': 'm06-2x',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Minnesota M06-HF Hybrid XC Functional': {
        'name': 'm06-hf',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'M06-L Meta-GGA XC Functional': {
        'name': 'm06-l',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Minnesota M08-HX Hybrid XC Functional': {
        'name': 'm08-hx',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'Minnesota M08-SO Hybrid XC Functional': {
        'name': 'm08-so',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'M11 Meta-GGA XC Functional': {
        'name': 'm11',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'M11-L Meta-GGA XC Functional': {
        'name': 'm11-l',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'MB3LYP-RC04 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mb3lyp-rc04',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MGGA_MS0 Meta-GGA XC Functional': {
        'name': 'mgga_ms0',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MGGA_MS1 Meta-GGA XC Functional': {
        'name': 'mgga_ms1',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MGGA_MS2 Meta-GGA XC Functional': {
        'name': 'mgga_ms2',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MGGA_MS2h Hybrid Meta-GGA XC Functional': {
        'name': 'mgga_ms2h',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MGGA_MVS Meta-GGA XC Functional': {
        'name': 'mgga_mvs',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MGGA_MV2h Hybrid Meta-GGA XC Functional': {
        'name': 'mgga_mvsh',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MN12-L Meta-GGA XC Functional': {
        'name': 'mn12-l',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'MN12-SX Meta-GGA Hybrid Screened Exchange-Correlation Functional': {
        'name': 'mn12-sx',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'MN15 Hybrid Meta-GGA Exchange-Correlation Functional': {
        'name': 'mn15',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'MN15-L Meta-GGA XC Functional': {
        'name': 'mn15-l',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MOHLYP GGA Exchange-Correlation Functional': {
        'name': 'mohlyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'MOHLYP2 GGA Exchange-Correlation Functional': {
        'name': 'mohlyp2',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPW1B95 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpw1b95',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'mPW1K Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpw1k',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPW1LYP Hybrid GGA Exchange-Correlation Functional': {
        'name': 'mpw1lyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPW1PBE Hybrid GGA Exchange-Correlation Functional': {
        'name': 'mpw1pbe',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPW1PW Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpw1pw',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'mPW3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpw3lyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPW3PW Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpw3pw',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPWB1K Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpwb1k',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'mPWLYP1M Hyb-GGA Exchange-Correlation Functional': {
        'name': 'mpwlyp1m',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPWLYP1W GGA Exchange-Correlation Functional': {
        'name': 'mpwlyp1w',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'mPWPW GGA Exchange-Correlation Functional': {
        'name': 'mpwpw',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'N12 nonseparable GGA Exchange-Correlation Functional': {
        'name': 'n12',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'N12-SX Hybrid nonseparable GGA Exchange-Correlation Functional': {
        'name': 'n12-sx',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'O3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'o3lyp',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'oblyp-e': {
        'name': 'oblyp-d',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE GGA Exchange, One-parameter Progressive correlation Functional': {
        'name': 'op-pbe',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'opbe-d': {
        'name': 'opbe-d',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'opwlyp-d': {
        'name': 'opwlyp-d',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'otpss-d': {
        'name': 'otpss-d',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE GGA Exchange-Correlation Functional': {
        'name': 'pbe',
        'dispersion': ['none', 'd3bj', 'd3mbj', 'nl'],
        'level': 'normal',
    },
    'PBE0 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'pbe0',
        'dispersion': ['none', 'd3bj', 'd3mbj', 'nl'],
        'level': 'normal',
    },
    'PBE0-13 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'pbe0-13',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE0-2 Double Hybrid Exchange-Correlation Functional': {
        'name': 'pbe0-2',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'pbe0-dh': {
        'name': 'pbe0-dh',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE1W GGA Exchange-Correlation Functional': {
        'name': 'pbe1w',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE50 Hybrid GGA Exchange-Correlation Functional': {
        'name': 'pbe50',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE Hybrid based 3C composite method with a small basis set, gCP and D3(BJ)': {  # noqa: E501
        'name': 'pbeh3c',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBELYP1W GGA Exchange-Correlation Functional': {
        'name': 'pbelyp1w',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PKZB Meta-GGA XC Functional': {
        'name': 'pkzb',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PTPSS SOS Double Hybrid XC Functional': {
        'name': 'ptpss',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'PW6B95 Hybrid Meta-GGA XC Functional': {
        'name': 'pw6b95',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'PW86B95 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'pw86b95',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PW86PBE GGA Exchange-Correlation Functional': {
        'name': 'pw86pbe',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PW91 GGA Exchange-Correlation Functional': {
        'name': 'pw91',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'PWB6K Hyb-GGA Exchange-Correlation Functional': {
        'name': 'pwb6k',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'PWPB95 SOS Double Hybrid XC Functional': {
        'name': 'pwpb95',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'revB3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'revb3lyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'revPBE GGA Exchange-Correlation Functional': {
        'name': 'revpbe',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'revPBE0 Hybrid GGA Exchange-Correlation Functional': {
        'name': 'revpbe0',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'revised TPSS Meta-GGA XC Functional': {
        'name': 'revtpss',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'revTPSSh Hyb-GGA Exchange-Correlation Functional': {
        'name': 'revtpssh',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'RPBE GGA Exchange-Correlation Functional': {
        'name': 'rpbe',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'SB98-1a Hyb-GGA Exchange-Correlation Functional': {
        'name': 'sb98-1a',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SB98-1b Hyb-GGA Exchange-Correlation Functional': {
        'name': 'sb98-1b',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SB98-1c Hyb-GGA Exchange-Correlation Functional': {
        'name': 'sb98-1c',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SB98-2a Hyb-GGA Exchange-Correlation Functional': {
        'name': 'sb98-2a',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SB98-2b Hyb-GGA Exchange-Correlation Functional': {
        'name': 'sb98-2b',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SB98-2c Hyb-GGA Exchange-Correlation Functional': {
        'name': 'sb98-2c',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SOGGA Exchange + PBE Correlation Functional': {
        'name': 'sogga',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SOGGA11 Exchange-Correlation Functional': {
        'name': 'sogga11',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'SOGGA11-X Hybrid Exchange-Correlation Functional': {
        'name': 'sogga11-x',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'SVWN3 (RPA) LSDA Functional': {
        'name': 'svwn',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TETER93 Exchange-Correlation Functional': {
        'name': 'teter93',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH-FC GGA Exchange-Correlation Functional': {
        'name': 'th-fc',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH-FCFO GGA Exchange-Correlation Functional': {
        'name': 'th-fcfo',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH-FCO GGA Exchange-Correlation Functional': {
        'name': 'th-fco',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH-FL GGA Exchange-Correlation Functional': {
        'name': 'th-fl',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH1 GGA Exchange-Correlation Functional': {
        'name': 'th1',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH2 GGA Exchange-Correlation Functional': {
        'name': 'th2',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH3 GGA Exchange-Correlation Functional': {
        'name': 'th3',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TH4 GGA Exchange-Correlation Functional': {
        'name': 'th4',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TPSS Meta-GGA XC Functional': {
        'name': 'tpss',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'TPSSh Hyb-GGA Exchange-Correlation Functional': {
        'name': 'tpssh',
        'dispersion': ['none', 'd3bj', 'nl'],
        'level': 'expert',
    },
    'TPSSLYP1W GGA Exchange-Correlation Functional': {
        'name': 'tpsslyp1w',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'TUNED-CAM-B3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'tuned-cam-b3lyp',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'VSXC Meta-GGA XC Functional': {
        'name': 'vsxc',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'VV10 GGA Exchange-Correlation Functional': {
        'name': 'vv10',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'wB97 GGA Exchange-Correlation Functional': {
        'name': 'wb97',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'wB97M-V Hyb-GGA Exchange-Correlation Functional': {
        'name': 'wb97m-v',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'wB97X Hyb-GGA Exchange-Correlation Functional': {
        'name': 'wb97x',
        'dispersion': ['none', 'd3bj', 'd'],
        'level': 'expert',
    },
    'wB97X-V Hyb-GGA Exchange-Correlation Functional': {
        'name': 'wb97x-v',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'PBE SR-XC Functional (HJS Model)': {
        'name': 'wpbe',
        'dispersion': ['none', 'd3bj', 'd3mbj'],
        'level': 'expert',
    },
    'PBE0 SR-XC Functional (HJS Model)': {
        'name': 'wpbe0',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'X1B95 Hyb-GGA Exchange-Correlation Functional': {
        'name': 'x1b95',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'X3LYP Hyb-GGA Exchange-Correlation Functional': {
        'name': 'x3lyp',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'XB1K Hyb-GGA Exchange-Correlation Functional': {
        'name': 'xb1k',
        'dispersion': ['none'],
        'level': 'expert',
    },
    'XLYP GGA Exchange-Correlation Functional': {
        'name': 'xlyp',
        'dispersion': ['none', 'd3bj'],
        'level': 'expert',
    },
    'ZLP GGA Exchange-Correlation Functional': {
        'name': 'zlp',
        'dispersion': ['none'],
        'level': 'expert',
    },
}  # yapf: disable

properties = {
    "-D ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "dispersion correction energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "Ha"
    },
    "32-POLE XXXXX": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXXXY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXXXZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXXYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXXYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXXZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXYYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXYYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXYZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XXZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XYYYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XYYYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XYYZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XYZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE XZZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE YYYYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE YYYYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE YYYZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE YYZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE YZZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "32-POLE ZZZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "5th order electrical multipole",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DFT FUNCTIONAL TOTAL ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "total energy of DFT functional",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DFT TOTAL ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "Total DFT energy including dispersion",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DFT VV10 ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "VV10 energy in DFT",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DFT XC ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "DFT exchange-correlation energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DIPOLE X": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical dipole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DIPOLE Y": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical dipole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DIPOLE Z": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical dipole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "DISPERSION CORRECTION ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "energy of the dispersion correction",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "ESP": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrostatic potential at nuclei",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XXXX": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XXXY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XXXZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XXYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XXYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XXZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XYYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XYYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XYZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE XZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE YYYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE YYYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE YYZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE YZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "HEXADECAPOLE ZZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical hexadecapole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "LOWDIN_CHARGES": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "atomic charges using Lowdin's method",
        "dimensionality": ["n_atoms"],
        "type": "float",
        "units": ""
    },
    "MAYER_INDICES": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "Mayer's bond indices",
        "dimensionality": ["triangular", "n_atoms", "n_atoms"],
        "type": "float",
        "units": ""
    },
    "MULLIKEN_CHARGES": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "atomic charges using Mulliken's method",
        "dimensionality": ["n_atoms"],
        "type": "float",
        "units": ""
    },
    "NUCLEAR REPULSION ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE XXX": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE XXY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE XXZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE XYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE XYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE XZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE YYY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE YYZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE YZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "OCTUPOLE ZZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical octupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "ONE-ELECTRON ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "the one-electron energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "PCM POLARIZATION ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "the PCM polarization energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "QUADRUPOLE XX": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical quadrupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "QUADRUPOLE XY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical quadrupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "QUADRUPOLE XZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical quadrupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "QUADRUPOLE YY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical quadrupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "QUADRUPOLE YZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical quadrupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "QUADRUPOLE ZZ": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical quadrupole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "SCF DIPOLE X": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical dipole moment from SCF",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "SCF DIPOLE Y": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical dipole moment from SCF",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "SCF DIPOLE Z": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "electrical dipole moment from SCF",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "SCF ITERATIONS": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "number of itereations in the SCF",
        "dimensionality": "scalar",
        "type": "integer",
        "units": ""
    },
    "TWO-ELECTRON ENERGY": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "the two-electron energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": ""
    },
    "WIBERG_LOWDIN_INDICES": {
        "calculation": [
            "energy",
            "optimization",
            "thermodynamics",
            "vibrations"
        ],
        "description": "the Wiber-Lowdin bond indices",
        "dimensionality": ["triangular", "n_atoms", "n_atoms"],
        "type": "float",
        "units": ""
    },
}  # yapf: disable
