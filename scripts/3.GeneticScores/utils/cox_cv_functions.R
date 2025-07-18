#' ----------------------------
#' Cox Proportional Hazards Cross-Validation Functions
#' 
#' Funciones reutilizables para entrenar y evaluar modelos CoxPH
#' con cross-validation en datos MDS.
#' 
#' Basado en la lógica de recomputeIPSSM.R
#' ----------------------------

library(survival)

#' Entrena modelo Cox en un fold específico
#' 
#' @param data Data frame con todas las variables
#' @param fold_id ID del fold a usar como entrenamiento (otros son test)
#' @param formula Formula del modelo Cox o vector de nombres de variables
#' @param time_var Nombre de la variable de tiempo (default: "LFS_YEARS")
#' @param status_var Nombre de la variable de estado (default: "LFS_STATUS")
#' @param fold_col Nombre de la columna de fold (default: "fold")
#' @return Lista con modelo entrenado y metadatos
train_cox_fold <- function(data, fold_id, formula, 
                          time_var = "LFS_YEARS", 
                          status_var = "LFS_STATUS",
                          fold_col = "fold") {
  
  # Filtrar datos de entrenamiento (excluir fold_id)
  train_data <- data[data[[fold_col]] != fold_id, ]
  
  # Crear formula si se pasa vector de variables
  if (is.character(formula)) {
    vars_str <- paste(formula, collapse = " + ")
    formula <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ ", vars_str))
  }
  
  # Entrenar modelo
  model <- coxph(formula, data = train_data)
  
  return(list(
    model = model,
    fold_id = fold_id,
    n_train = nrow(train_data),
    formula = formula
  ))
}

#' Evalúa modelo Cox en datos de test
#' 
#' @param trained_model Lista resultado de train_cox_fold()
#' @param test_data Data frame con datos de test
#' @param time_var Nombre de la variable de tiempo
#' @param status_var Nombre de la variable de estado
#' @return Lista con métricas de evaluación
evaluate_cox_fold <- function(trained_model, test_data,
                             time_var = "LFS_YEARS",
                             status_var = "LFS_STATUS") {
  
  model <- trained_model$model
  
  # Generar predicciones (log hazard ratio)
  predictions <- predict(model, newdata = test_data, type = "lp")
  
  # Calcular C-index
  surv_obj <- Surv(test_data[[time_var]], test_data[[status_var]])
  predict_conc <- -predictions  # Invertir para concordancia
  c_index <- concordance(surv_obj ~ predict_conc)$concordance

  return(list(
    fold_id = trained_model$fold_id,
    c_index = c_index,
    predictions = predictions,
    n_test = nrow(test_data),
    n_events = sum(test_data[[status_var]])
  ))
}

#' Ejecuta cross-validation completo para modelo Cox
#' 
#' @param data Data frame con todas las variables incluyendo folds
#' @param formula Formula del modelo o vector de nombres de variables
#' @param time_var Nombre de la variable de tiempo
#' @param status_var Nombre de la variable de estado
#' @param fold_col Nombre de la columna de fold
#' @return Lista con resultados agregados de CV
run_cox_cv <- function(data, formula,
                      time_var = "LFS_YEARS",
                      status_var = "LFS_STATUS", 
                      fold_col = "fold") {
  
  # Obtener folds únicos
  unique_folds <- unique(data[[fold_col]])
  unique_folds <- unique_folds[!is.na(unique_folds)]
  
  results <- list()
  all_predictions <- c()
  all_truth <- list()
  
  # Iterar sobre cada fold
  for (fold_id in unique_folds) {
    
    # Entrenar modelo
    trained_model <- train_cox_fold(data, fold_id, formula, time_var, status_var, fold_col)
    
    # Datos de test (solo el fold actual)
    test_data <- data[data[[fold_col]] == fold_id, ]
    
    # Evaluar modelo
    fold_result <- evaluate_cox_fold(trained_model, test_data, time_var, status_var)
    
    # Guardar resultados
    results[[paste0("fold_", fold_id)]] <- fold_result
    
    # Acumular predicciones para C-index global
    all_predictions <- c(all_predictions, fold_result$predictions)
    all_truth$time <- c(all_truth$time, test_data[[time_var]])
    all_truth$status <- c(all_truth$status, test_data[[status_var]])
  }
  
  # Calcular C-index global
  global_surv <- Surv(all_truth$time, all_truth$status)
  all_predic_concord <- -all_predictions  # Invertir para concordancia
  global_c_index <- concordance(global_surv ~ all_predic_concord)$concordance

  # Calcular estadísticas agregadas
  fold_c_indices <- sapply(results, function(x) x$c_index)
  
  return(list(
    fold_results = results,
    global_c_index = global_c_index,
    mean_c_index = mean(fold_c_indices, na.rm = TRUE),
    sd_c_index = sd(fold_c_indices, na.rm = TRUE),
    n_folds = length(unique_folds),
    formula = formula
  ))
}
