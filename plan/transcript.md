# Prompts used in this development

## Prompt 1

You are working in the lr_reduction project, in the bvacaliuc/new_workflow_ui_plan branch. Examine the git log on this branch and consider the files added relative to the next branch. This branch has implemented a new set of functions that use reduce_from_template(). There is an example in example_nr_from_template.py that demonstrates how to setup and call this function. All the needed datasets are provided in the system at /SNS/REF_L/IPTS-30101/ and /SNS/REF_L/shared/. Review the code in launcher/launcher.py. Your goal is to develop a plan to extend launcher.py so that it has the following feature: A new tab that exposes all the possible configuration options available in reduce_from_template(). This tab has a 'REDUCE' button that when pressed will invoke reduce_from_template() gathering the options in the tab and passing them to reduce_from_template(). Every path selector uses a folder browser with the expected file extension filter. Add a progress bar to the tab to indicate how far the reduction has progressed and when it is complete. Make the necessary changes to the code files to support this. Make sure the output folder is a path selector. Add the option to save plots, and if selected, ensure the plots are written directly to files in the output_dir. Add an option to interact with the plots; if the user chooses not to interact, do not show the plot. Please create a detailed plan for this work. I would like to review the development plan before you implement it.

The plan that Claude developed from this prompt is [add-template-reduction-tab-to-launcher.md](add-template-reduction-tab-to-launcher.md)

### Prompt 1.1

When starting the launcher via 'pixi run python launcher/launcher.py', there is a missing dependency. Please review and update. Please develop a set of comprehensive UI tests with offline rendering now that the launcher is being integrated into the backend. Also review the pre-existing errors in test_web_report.py and explain to me their nature. Give me some options on what the fixes would require.

### Prompt 1.2

Please write your analysis above to a plan that can be acted upon in a future session.

The plan that Claude developed is [fix-test-web-report-errors.md](fix-test-web-report-errors.md).

Since Option D: leave as-is will not break CI, we can defer engaging this to a separate step.

### Prompt 1.3

When running 'pixi run python launcher/launcher.py' and pressing REDUCE, gives a fault. Please investigate that fault and consider as many other possible failure modes that you can due to UI elements. It is imperative that the UI have meaningful data entry validation and be robust. The test system should be expanded to cover such UI elements. Please add an option to log output to stdout for aid in diagnosing such errors.

### Prompt 1.4

I get this error when pressing REDUCE:
```
--- Template Reduction ---
  Run numbers: [211029]
  Template: /SNS/REF_L/IPTS-30101/shared/autoreduce/REF_L_211031_auto_template.xml
  Experiment: IPTS-30101
  Data path: /SNS/REF_L/IPTS-30101/nexus
  Template save: /tmp/IPTS-30101
  Save plots: True  Plot dir: /tmp/IPTS-30101
  Override params: {'method': 'meantheta', 'Normalize': True, 'AutoScale': True, 'useCalcTheta': False, 'plotQ4': False, 'DetResFn': 'rectangular', 'peak_pad': 1, 'peak_type': 'supergauss', 'use_emission_time': True, 'Spath': '/tmp/IPTS-30101', 'DBpath': '/SNS/REF_L/shared/Cd_DB_processing/DBs', 'DBname': ['A1_div5_Cd_DB.dat']}
--------------------------
[Worker] Starting run 211029 (1/1)
[Worker] Run 211029: Reading NEXUS metadata (0%)
[Worker] ERROR:
Traceback (most recent call last):
  File "/home/bvacaliuc/Projects/Claude/lr_reduction/launcher/apps/reduction_worker.py", line 60, in run
    result = reduce_from_template(
             ^^^^^^^^^^^^^^^^^^^^^
  File "/home/bvacaliuc/Projects/Claude/lr_reduction/src/lr_reduction/new_reduction_from_template.py", line 56, in reduce_from_template
    fname = datapath / f"REF_L_{runno}.nxs.h5"
            ~~~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~
TypeError: unsupported operand type(s) for /: 'str' and 'str'

REDUCTION ERROR:
Traceback (most recent call last):
  File "/home/bvacaliuc/Projects/Claude/lr_reduction/launcher/apps/reduction_worker.py", line 60, in run
    result = reduce_from_template(
             ^^^^^^^^^^^^^^^^^^^^^
  File "/home/bvacaliuc/Projects/Claude/lr_reduction/src/lr_reduction/new_reduction_from_template.py", line 56, in reduce_from_template
    fname = datapath / f"REF_L_{runno}.nxs.h5"
            ~~~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~
TypeError: unsupported operand type(s) for /: 'str' and 'str'
```

### Prompt 1.5

I ran the REDUCE successfully. Please check on these RuntimeWarnings that were emitted and improve handling of such conditions. This is the output I got:
```
--- Template Reduction ---
  Run numbers: [211029, 211030, 211031]
  Template: /SNS/REF_L/IPTS-30101/shared/autoreduce/REF_L_211031_auto_template.xml
  Experiment: IPTS-30101
  Data path: /SNS/REF_L/IPTS-30101/nexus
  Template save: /tmp/IPTS-30101
  Save plots: True  Plot dir: /tmp/IPTS-30101
  Override params: {'method': 'meantheta', 'Normalize': True, 'AutoScale': True, 'useCalcTheta': False, 'plotQ4': False, 'DetResFn': 'rectangular', 'peak_pad': 1, 'peak_type': 'supergauss', 'use_emission_time': True, 'Spath': '/tmp/IPTS-30101', 'DBpath': '/SNS/REF_L/shared/Cd_DB_processing/DBs', 'DBname': ['A1_div5_Cd_DB.dat']}
--------------------------
[Worker] Starting run 211029 (1/3)
[Worker] Run 211029: Reading NEXUS metadata (0%)
[Worker] Run 211029: Configuring reduction (6%)
[Worker] Run 211029: Running reduction (13%)
Older run doesn't include coordinates PV
Binary data computed successfully for run 211029
Theta, TTHD = 0.515, 1.0
Calculated theta: 0.614, dTheta: 0.098
meantheta
Normalization factor: 0.089
Saved combined result to /tmp/IPTS-30101/REFL_211029_1_211029_partial.dat
[Worker] Run 211029: Assembling results (20%)
Files found: 0
Data loaded: 1
Saved combined result to /tmp/IPTS-30101/REFL_211029_combined_data.dat
[Worker] Run 211029: Saving template (26%)
Reading file /SNS/REF_L/IPTS-30101/shared/autoreduce/REF_L_211031_auto_template.xml
Saving to file REFL_211029_template_new.xml
[Worker] Starting run 211030 (2/3)
[Worker] Run 211030: Reading NEXUS metadata (33%)
[Worker] Run 211030: Configuring reduction (40%)
[Worker] Run 211030: Running reduction (46%)
Older run doesn't include coordinates PV
Binary data computed successfully for run 211030
Theta, TTHD = 1.316, 2.601
Calculated theta: 1.419, dTheta: 0.103
meantheta
/home/bvacaliuc/Projects/Claude/lr_reduction/.pixi/envs/default/lib/python3.11/site-packages/numpy/_core/fromnumeric.py:3904: RuntimeWarning: Mean of empty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
/home/bvacaliuc/Projects/Claude/lr_reduction/.pixi/envs/default/lib/python3.11/site-packages/numpy/_core/_methods.py:147: RuntimeWarning: invalid value encountered in scalar divide
  ret = ret.dtype.type(ret / rcount)
Normalization factor: nan
Saved combined result to /tmp/IPTS-30101/REFL_211029_2_211030_partial.dat
[Worker] Run 211030: Assembling results (53%)
Files found: 0
Data loaded: 2
/home/bvacaliuc/Projects/Claude/lr_reduction/src/lr_reduction/nr_tools.py:76: RuntimeWarning: invalid value encountered in scalar divide
  mean=np.sum(v*w)/np.sum(w)
/home/bvacaliuc/Projects/Claude/lr_reduction/src/lr_reduction/nr_tools.py:84: RuntimeWarning: invalid value encountered in scalar divide
  mean = np.sum(v * w) / np.sum(w)
/home/bvacaliuc/Projects/Claude/lr_reduction/src/lr_reduction/nr_tools.py:85: RuntimeWarning: divide by zero encountered in scalar divide
  sigma_mean = np.sqrt(1 / np.sum(w))
Scaling factor: nan
Saved combined result to /tmp/IPTS-30101/REFL_211029_combined_data.dat
[Worker] Run 211030: Saving template (60%)
Reading file /tmp/IPTS-30101/REFL_211029_template_new.xml
Using prior template data
Saving to file REFL_211029_template_new.xml
[Worker] Starting run 211031 (3/3)
[Worker] Run 211031: Reading NEXUS metadata (66%)
[Worker] Run 211031: Configuring reduction (73%)
[Worker] Run 211031: Running reduction (80%)
Older run doesn't include coordinates PV
Binary data computed successfully for run 211031
Theta, TTHD = 3.518, 7.001
Calculated theta: 3.624, dTheta: 0.106
meantheta
Normalization factor: nan
Saved combined result to /tmp/IPTS-30101/REFL_211029_3_211031_partial.dat
[Worker] Run 211031: Assembling results (86%)
Files found: 0
Data loaded: 3
Scaling factor: nan
Scaling factor: nan
Saved combined result to /tmp/IPTS-30101/REFL_211029_combined_data.dat
[Worker] Run 211031: Saving template (93%)
Reading file /tmp/IPTS-30101/REFL_211029_template_new.xml
Using prior template data
Saving to file REFL_211029_template_new.xml
[Worker] All 3 run(s) completed successfully
```

## Prompt 2

You are working in the lr_reduction project, in the bvacaliuc/new_workflow_ui_plan branch. There is a deficiency in the way that pixi is handled when the work is shared among multiple users. If one user runs pixi install, other users are prevented from using pixi run the folder due to permissions problems on the pixi state. Examine the nsd-app-wrap project for a key step with respect to pixi deployments that has the necessary syntax the developers have worked out to solve this problem on analysis.sns.gov and the deployments made to /usr/local/pixi/. The goal of this task is to understand the needed adjustments to managing the pixi environment properly to avoid the need of writing to shared files in the .pixi/ of the project when doing 'pixi run {command} ...'. I would like to see your plan before you implement it.

The explanation and plan that Claude came up with is [fix-shared-user-pixi-permission.md](fix-shared-user-pixi-permission.md)

### Prompt 2.1

Regarding steps 2 thru 4, these steps are done when the project is deployed. My concern is when multiple developers are working in a shared space prior to deployment. I would like the instructions modified to detail the way that users can interact collaboratively with a common shared clone of the project. Is it possible? Feasible? Recommended? Help me to understand my options here.

