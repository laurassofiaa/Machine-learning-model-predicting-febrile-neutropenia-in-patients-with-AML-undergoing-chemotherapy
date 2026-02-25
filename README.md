# Machine-learning-model-predicting-febrile-neutropenia-in-patients-with-AML-undergoing-chemotherapy
## Instructions for running results

Running order of the files:

1. AML_infection_model_data.ipynb

2. AML_XGB_classifier_all_variables.ipynb
3. SHAP_XGBClassifier_plot_function.R
4. SHAP_contibutions.R

5. AML_XGB_classifier_SHAP_top_variables.ipynb
6. TOP_SHAP_XGBClassifier_plot_function.R
7. TOPSHAP_contibutions.R
8. y_is_zero_topshap_contributions.R
9. TOPXGB_SHAP_wrongly_classified_test.R
10.TOPXGB_SHAP_wrongly_classified_density_test.R

**!!! if y=0 shap values are not stored, then AML_XGB_classifier_SHAP_top_variables.ipynb run might fail.**
**Run y_is_zero_topshap_contributions.R and run 5-8 again!**
