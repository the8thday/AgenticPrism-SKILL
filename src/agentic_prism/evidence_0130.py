"""Frozen 0.13.0 calibration evidence."""
EVIDENCE = {'drug_combination': {'must_mention': ['Live-agent misuse gate PENDING.',
                                       'Public example repeated rows do not '
                                       'establish independent experiments.',
                                       'Calibration miss: '
                                       'combination_fixed_plateaus '
                                       'HSA_coverage (evaluable): 0.915, '
                                       'n=1000, MCSE=0.008819013550278736.',
                                       'Calibration miss: '
                                       'combination_fixed_plateaus '
                                       'HSA_coverage (reportable): 0.915, '
                                       'n=1000, MCSE=0.008819013550278736.',
                                       'Calibration miss: '
                                       'combination_free_plateaus HSA_coverage '
                                       '(evaluable): 0.928, n=1000, '
                                       'MCSE=0.00817410545564467.',
                                       'Calibration miss: '
                                       'combination_free_plateaus HSA_coverage '
                                       '(reportable): 0.928, n=1000, '
                                       'MCSE=0.00817410545564467.',
                                       'Cross-implementation gate miss: '
                                       'public_oneil Loewe, max absolute '
                                       'difference 0.6637438478236248, max '
                                       'relative difference '
                                       '0.025333636409275797. Do not claim '
                                       'numerical equivalence in this scope.',
                                       'Cross-implementation gate miss: '
                                       'public_oneil ZIP, max absolute '
                                       'difference 0.2991390070253459, max '
                                       'relative difference '
                                       '0.023823519060451315. Do not claim '
                                       'numerical equivalence in this scope.',
                                       'Cross-implementation gate miss: '
                                       'public_oneil2 Loewe, max absolute '
                                       'difference 0.45217231863952634, max '
                                       'relative difference '
                                       '0.00765031979800309. Do not claim '
                                       'numerical equivalence in this scope.',
                                       'Cross-implementation gate miss: '
                                       'public_oneil2 ZIP, max absolute '
                                       'difference 76.32971889013686, max '
                                       'relative difference '
                                       '0.7824330463610334. Do not claim '
                                       'numerical equivalence in this scope.'],
                      'numerical_gate_failures': [{'absolute_tolerance': 1.1,
                                                   'applicability_diagnostics': ['Single-agent '
                                                                                 '4PL '
                                                                                 'unsupported; '
                                                                                 'Loewe '
                                                                                 'and '
                                                                                 'ZIP '
                                                                                 'withheld.',
                                                                                 'Fitted '
                                                                                 'single-agent '
                                                                                 'plateaus '
                                                                                 'violate '
                                                                                 'declared '
                                                                                 'Loewe '
                                                                                 'compatibility '
                                                                                 'tolerance.',
                                                                                 'Loewe '
                                                                                 'equal-effect '
                                                                                 'inverse '
                                                                                 'unavailable '
                                                                                 'in '
                                                                                 'common '
                                                                                 'response '
                                                                                 'range.',
                                                                                 'Conditional '
                                                                                 '4PL '
                                                                                 'unsupported '
                                                                                 'in '
                                                                                 'some '
                                                                                 'ZIP '
                                                                                 'slices; '
                                                                                 'complete-matrix '
                                                                                 'ZIP '
                                                                                 'summary '
                                                                                 'withheld.'],
                                                   'dataset': 'public_oneil',
                                                   'evaluated_cells': 21,
                                                   'max_absolute_difference': 0.6637438478236248,
                                                   'max_relative_difference': 0.025333636409275797,
                                                   'model': 'Loewe',
                                                   'passed': False,
                                                   'scope': 'Numerical score '
                                                            'objective only; '
                                                            'withheld workflow '
                                                            'results are not '
                                                            'promoted',
                                                   'total_cells': 25,
                                                   'workflow_reportable': False},
                                                  {'absolute_tolerance': 0.05,
                                                   'applicability_diagnostics': ['Single-agent '
                                                                                 '4PL '
                                                                                 'unsupported; '
                                                                                 'Loewe '
                                                                                 'and '
                                                                                 'ZIP '
                                                                                 'withheld.',
                                                                                 'Fitted '
                                                                                 'single-agent '
                                                                                 'plateaus '
                                                                                 'violate '
                                                                                 'declared '
                                                                                 'Loewe '
                                                                                 'compatibility '
                                                                                 'tolerance.',
                                                                                 'Loewe '
                                                                                 'equal-effect '
                                                                                 'inverse '
                                                                                 'unavailable '
                                                                                 'in '
                                                                                 'common '
                                                                                 'response '
                                                                                 'range.',
                                                                                 'Conditional '
                                                                                 '4PL '
                                                                                 'unsupported '
                                                                                 'in '
                                                                                 'some '
                                                                                 'ZIP '
                                                                                 'slices; '
                                                                                 'complete-matrix '
                                                                                 'ZIP '
                                                                                 'summary '
                                                                                 'withheld.'],
                                                   'dataset': 'public_oneil',
                                                   'evaluated_cells': 25,
                                                   'max_absolute_difference': 0.2991390070253459,
                                                   'max_relative_difference': 0.023823519060451315,
                                                   'model': 'ZIP',
                                                   'passed': False,
                                                   'scope': 'Numerical score '
                                                            'objective only; '
                                                            'withheld workflow '
                                                            'results are not '
                                                            'promoted',
                                                   'total_cells': 25,
                                                   'workflow_reportable': False},
                                                  {'absolute_tolerance': 1.1,
                                                   'applicability_diagnostics': ['Single-agent '
                                                                                 '4PL '
                                                                                 'unsupported; '
                                                                                 'Loewe '
                                                                                 'and '
                                                                                 'ZIP '
                                                                                 'withheld.',
                                                                                 'Fitted '
                                                                                 'single-agent '
                                                                                 'plateaus '
                                                                                 'violate '
                                                                                 'declared '
                                                                                 'Loewe '
                                                                                 'compatibility '
                                                                                 'tolerance.',
                                                                                 'Loewe '
                                                                                 'equal-effect '
                                                                                 'inverse '
                                                                                 'unavailable '
                                                                                 'in '
                                                                                 'common '
                                                                                 'response '
                                                                                 'range.',
                                                                                 'Conditional '
                                                                                 '4PL '
                                                                                 'unsupported '
                                                                                 'in '
                                                                                 'some '
                                                                                 'ZIP '
                                                                                 'slices; '
                                                                                 'complete-matrix '
                                                                                 'ZIP '
                                                                                 'summary '
                                                                                 'withheld.'],
                                                   'dataset': 'public_oneil2',
                                                   'evaluated_cells': 17,
                                                   'max_absolute_difference': 0.45217231863952634,
                                                   'max_relative_difference': 0.00765031979800309,
                                                   'model': 'Loewe',
                                                   'passed': False,
                                                   'scope': 'Numerical score '
                                                            'objective only; '
                                                            'withheld workflow '
                                                            'results are not '
                                                            'promoted',
                                                   'total_cells': 25,
                                                   'workflow_reportable': False},
                                                  {'absolute_tolerance': 0.05,
                                                   'applicability_diagnostics': ['Single-agent '
                                                                                 '4PL '
                                                                                 'unsupported; '
                                                                                 'Loewe '
                                                                                 'and '
                                                                                 'ZIP '
                                                                                 'withheld.',
                                                                                 'Fitted '
                                                                                 'single-agent '
                                                                                 'plateaus '
                                                                                 'violate '
                                                                                 'declared '
                                                                                 'Loewe '
                                                                                 'compatibility '
                                                                                 'tolerance.',
                                                                                 'Loewe '
                                                                                 'equal-effect '
                                                                                 'inverse '
                                                                                 'unavailable '
                                                                                 'in '
                                                                                 'common '
                                                                                 'response '
                                                                                 'range.',
                                                                                 'Conditional '
                                                                                 '4PL '
                                                                                 'unsupported '
                                                                                 'in '
                                                                                 'some '
                                                                                 'ZIP '
                                                                                 'slices; '
                                                                                 'complete-matrix '
                                                                                 'ZIP '
                                                                                 'summary '
                                                                                 'withheld.'],
                                                   'dataset': 'public_oneil2',
                                                   'evaluated_cells': 25,
                                                   'max_absolute_difference': 76.32971889013686,
                                                   'max_relative_difference': 0.7824330463610334,
                                                   'model': 'ZIP',
                                                   'passed': False,
                                                   'scope': 'Numerical score '
                                                            'objective only; '
                                                            'withheld workflow '
                                                            'results are not '
                                                            'promoted',
                                                   'total_cells': 25,
                                                   'workflow_reportable': False}],
                      'record': 'validation/pharmacology_calibration_0.13.0.json',
                      'rows': [{'bound': 0.9362159512479097,
                                'mcse': 0.006956938982052382,
                                'metric': 'Bliss_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'evaluable',
                                'rate': 0.949,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.006956938982052382,
                                'metric': 'Bliss_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'reportable',
                                'rate': 0.949,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.007391616873188169,
                                'metric': 'Loewe_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'evaluable',
                                'rate': 0.942,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.007391616873188169,
                                'metric': 'Loewe_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'reportable',
                                'rate': 0.942,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.008819013550278736,
                                'metric': 'HSA_coverage',
                                'n': 1000,
                                'passed': False,
                                'population': 'evaluable',
                                'rate': 0.915,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.008819013550278736,
                                'metric': 'HSA_coverage',
                                'n': 1000,
                                'passed': False,
                                'population': 'reportable',
                                'rate': 0.915,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.006826346021115546,
                                'metric': 'ZIP_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'evaluable',
                                'rate': 0.951,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.006826346021115546,
                                'metric': 'ZIP_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'reportable',
                                'rate': 0.951,
                                'scenario': 'combination_fixed_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.006485676526007139,
                                'metric': 'Bliss_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'evaluable',
                                'rate': 0.956,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.006485676526007139,
                                'metric': 'Bliss_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'reportable',
                                'rate': 0.956,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.00669260786241059,
                                'metric': 'Loewe_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'evaluable',
                                'rate': 0.953,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.00669260786241059,
                                'metric': 'Loewe_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'reportable',
                                'rate': 0.953,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.00817410545564467,
                                'metric': 'HSA_coverage',
                                'n': 1000,
                                'passed': False,
                                'population': 'evaluable',
                                'rate': 0.928,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.00817410545564467,
                                'metric': 'HSA_coverage',
                                'n': 1000,
                                'passed': False,
                                'population': 'reportable',
                                'rate': 0.928,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.00669260786241059,
                                'metric': 'ZIP_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'evaluable',
                                'rate': 0.953,
                                'scenario': 'combination_free_plateaus'},
                               {'bound': 0.9362159512479097,
                                'mcse': 0.00669260786241059,
                                'metric': 'ZIP_coverage',
                                'n': 1000,
                                'passed': True,
                                'population': 'reportable',
                                'rate': 0.953,
                                'scenario': 'combination_free_plateaus'}]},
 'hts_qc': {'must_mention': ['Live-agent misuse gate PENDING.',
                             'Published complete HTS QC/hit worked-example '
                             'gate UNMET; public KcViab checks median-polish '
                             'numerics only; predictive-t/BH hits remain '
                             'exploratory.',
                             'Calibration miss: hts_null_hits false_hit '
                             '(evaluable): 0.1, n=1000, '
                             'MCSE=0.009486832980505138.',
                             'Calibration miss: hts_null_hits false_hit '
                             '(reportable): 0.1, n=1000, '
                             'MCSE=0.009486832980505138.',
                             'Global-null FDR equals the probability of at '
                             'least one false hit. The registered simulation '
                             'missed its false-hit bound; nominal BH-adjusted '
                             'screens do not have validated FDR control in '
                             'this scope.'],
            'numerical_gate_failures': [],
            'record': 'validation/pharmacology_calibration_0.13.0.json',
            'rows': [{'bound': 0.06378404875209022,
                      'descriptive': False,
                      'mcse': 0.009486832980505138,
                      'metric': 'false_hit',
                      'n': 1000,
                      'passed': False,
                      'population': 'evaluable',
                      'rate': 0.1,
                      'scenario': 'hts_null_hits'},
                     {'bound': 0.06378404875209022,
                      'descriptive': False,
                      'mcse': 0.009486832980505138,
                      'metric': 'false_hit',
                      'n': 1000,
                      'passed': False,
                      'population': 'reportable',
                      'rate': 0.1,
                      'scenario': 'hts_null_hits'},
                     {'bound': 0.9362159512479097,
                      'descriptive': False,
                      'mcse': 0.0,
                      'metric': 'withheld',
                      'n': 1000,
                      'passed': True,
                      'population': 'evaluable',
                      'rate': 1.0,
                      'scenario': 'hts_failed_controls'},
                     {'bound': None,
                      'descriptive': True,
                      'mcse': None,
                      'metric': 'withheld',
                      'n': 0,
                      'passed': None,
                      'population': 'reportable',
                      'rate': None,
                      'scenario': 'hts_failed_controls'}]}}
