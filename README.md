# iterative_bspline_model_fitter_for_quasar_continuum
The iterative B-spline model fitter for quasar continuum is a Python code adapted from the SDSS BOSS reduction pipeline originally written in IDL. The code leverages `pyidl` and `pyidlutils` to translate the IDL code into Python, with necessary adjustments to address differences in how numpy and IDL handle floating-point numbers. Specifically, modifications were made to eliminate any breakpoints with a difference of only one index.

## Logic Behind the Code

The code operates on a straightforward logic to determine breakpoints for fitting a B-spline model to the quasar continuum. Here’s a detailed explanation of the process:

1. **Weight Array Creation**: 
   - The flux values are used to create a weight array, which exaggerates regions where emission lines are present. This helps in assigning higher importance to these regions when determining breakpoints.
   - For regions with absorption lines, the weights are comparatively lower, ensuring breakpoints are more spaced out.

2. **Normalization**:
   - The weight array is normalized to ensure it sums up to 1, making the subsequent steps easier to manage and interpret.

3. **Cumulative Sum Array**:
   - A cumulative sum of the normalized weight array is calculated. This cumulative sum array helps in identifying where significant changes in the flux occur, which is critical for placing breakpoints.

4. **Integer Conversion and Index Finding**:
   - The cumulative sum array is converted to an integer array. This conversion effectively bins the continuous cumulative values into discrete steps.
   - The indices where the integer array changes its value are identified. These indices represent the breakpoints for the B-spline fitting. Breakpoints are closely spaced in regions with emission lines and more spaced out in regions with absorption lines.

5. **Refinement of Breakpoints**:
   - To ensure the breakpoints are well-placed, any breakpoints that are too close to each other (differing by only one index) are removed. This refinement step ensures a more stable and meaningful fit for the continuum.

6. **Model Spectrum**:
   - Using the identified breakpoints, a B-spline model is fitted to the quasar continuum. This model spectrum aims to accurately represent the continuum, distinguishing between regions of emission and absorption.
