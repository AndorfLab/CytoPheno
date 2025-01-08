# Input choices in the Shiny UI

## Tab 1
# Vector of choices for Median Difference Equation parameters
input_choices_marker_diff_eq <- c("default_marker_diff_eq", "choose_marker_diff_eq")
names(input_choices_marker_diff_eq) <- c("Yes", "No")

# Vector of choices for transform data
input_choices_transform <- c("transform_no", "transform_yes")
names(input_choices_transform) <- c("No", "Yes")

# Vector of choices for cell type reference
input_choices_ref_1 <- c("internal_ref_1", "upload_ref_1")
names(input_choices_ref_1) <- c("Default references", "File upload")

# Vector of choices for cell type ontology
input_choices_onto_1 <- c("internal_CL_1", "internal_pCL_1")
names(input_choices_onto_1) <- c("Cell ontology", "Provisional cell ontology")

# Vector of choices for "high"
input_choices_high_1 <- c("high_high_1","high_positive_1")
names(input_choices_high_1) <- c("Only include 'high' matches","Include 'positive' marker signs when searching 'high'")

# Vector of choices for "low"
input_choices_low_1 <- c("low_low_1","low_positive_1","low_negative_1")
names(input_choices_low_1) <- c("Only include 'low' matches","Include 'positive' marker signs when searching 'low'", "Include 'negative' marker signs when searching 'low'")

#### Tab 2
# Vector of choices for first button (file type)
input_choices <- c("type_input", "csv_upload")
names(input_choices) <- c("Text input", "File upload")

# Vector of choices for cell type reference
input_choices_ref_2 <- c("internal_ref_2", "upload_ref_2")
names(input_choices_ref_2) <- c("Default references", "File upload")

# Vector of choices for cell type ontology
input_choices_onto_2 <- c("internal_CL_2", "internal_pCL_2")
names(input_choices_onto_2) <- c("Cell ontology", "Provisional cell ontology")

# Vector of choices for "high"
input_choices_high_upload <- c("high_high_upload","high_positive_upload")
names(input_choices_high_upload) <- c("Only include 'high' matches","Include 'positive' marker signs when searching 'high'")

input_choices_high_2 <- c("high_high_2","high_positive_2")
names(input_choices_high_2) <- c("Only include 'high' matches","Include 'positive' marker signs when searching 'high'")

# Vector of choices for "low"
input_choices_low_upload <- c("low_low_upload","low_positive_upload","low_negative_upload")
names(input_choices_low_upload) <- c("Only include 'low' matches","Include 'positive' marker signs when searching 'low'", "Include 'negative' marker signs when searching 'low'")

input_choices_low_2 <- c("low_low_2","low_positive_2","low_negative_2")
names(input_choices_low_2) <- c("Only include 'low' matches","Include 'positive' marker signs when searching 'low'", "Include 'negative' marker signs when searching 'low'")
