#' ----------------------------
#' Tests para Cox CV Functions
#' 
#' Tests unitarios para verificar que las funciones de cross-validation
#' de Cox funcionan correctamente con datos sintéticos.
#' ----------------------------

# Cargar funciones a testear
source("../utils/cox_cv_functions.R")

# Función para crear datos sintéticos pequeños
create_synthetic_data <- function(n = 50, n_folds = 3) {
  set.seed(123)
  
  data.frame(
    ID = 1:n,
    AGE = rnorm(n, 65, 10),
    SEX = sample(c("M", "F"), n, replace = TRUE),
    HB = rnorm(n, 10, 2),
    PLT = rnorm(n, 200, 50),
    TP53 = rbinom(n, 1, 0.3),
    ASXL1 = rbinom(n, 1, 0.2),
    LFS_YEARS = rexp(n, 0.1),
    LFS_STATUS = rbinom(n, 1, 0.4),
    fold = sample(1:n_folds, n, replace = TRUE)
  )
}

# Test 1: Verificar que train_cox_fold no falla
test_train_cox_fold_basic <- function() {
  cat("Test 1: train_cox_fold - funcionamiento básico\n")
  
  # Datos sintéticos
  data <- create_synthetic_data()
  
  # Variables a usar
  vars <- c("AGE", "SEX", "HB", "TP53")
  
  # Entrenar modelo en fold 1
  result <- train_cox_fold(data, fold_id = 1, formula = vars)
  
  # Verificaciones básicas
  stopifnot("model" %in% names(result))
  stopifnot("fold_id" %in% names(result))
  stopifnot("n_train" %in% names(result))
  stopifnot("formula" %in% names(result))
  
  # Verificar que el modelo es de clase coxph
  stopifnot(inherits(result$model, "coxph"))
  
  # Verificar que fold_id es correcto
  stopifnot(result$fold_id == 1)
  
  # Verificar que n_train es razonable
  stopifnot(result$n_train > 0)
  stopifnot(result$n_train < nrow(data))  # Debe ser menor que total
  
  cat("  ✓ Función ejecuta sin errores\n")
  cat("  ✓ Formato de output correcto\n")
  cat("  ✓ Modelo Cox creado correctamente\n")
  cat("  ✓ Metadatos incluidos\n")
}

# Test 2: Verificar que evaluate_cox_fold no falla
test_evaluate_cox_fold_basic <- function() {
  cat("\nTest 2: evaluate_cox_fold - funcionamiento básico\n")
  
  # Datos sintéticos
  data <- create_synthetic_data()
  vars <- c("AGE", "HB", "TP53")
  
  # Entrenar modelo
  trained_model <- train_cox_fold(data, fold_id = 1, formula = vars)
  
  # Datos de test (fold 1)
  test_data <- data[data$fold == 1, ]
  
  # Evaluar modelo
  result <- evaluate_cox_fold(trained_model, test_data)
  
  # Verificaciones básicas
  stopifnot("fold_id" %in% names(result))
  stopifnot("c_index" %in% names(result))
  stopifnot("predictions" %in% names(result))
  stopifnot("n_test" %in% names(result))
  stopifnot("n_events" %in% names(result))
  
  # Verificar tipos de datos
  stopifnot(is.numeric(result$c_index))
  stopifnot(is.numeric(result$predictions))
  stopifnot(is.numeric(result$n_test))
  stopifnot(is.numeric(result$n_events))
  
  # Verificar rangos razonables
  stopifnot(result$c_index >= 0 && result$c_index <= 1)
  stopifnot(result$n_test == nrow(test_data))
  stopifnot(result$n_events <= result$n_test)
  stopifnot(length(result$predictions) == result$n_test)
  
  cat("  ✓ Función ejecuta sin errores\n")
  cat("  ✓ Formato de output correcto\n")
  cat("  ✓ C-index en rango válido [0,1]\n")
  cat("  ✓ Dimensiones consistentes\n")
}

# Test 3: Verificar que run_cox_cv completo no falla
test_run_cox_cv_basic <- function() {
  cat("\nTest 3: run_cox_cv - funcionamiento completo\n")
  
  # Datos sintéticos
  data <- create_synthetic_data(n = 60, n_folds = 4)
  vars <- c("AGE", "HB", "TP53", "ASXL1")
  
  # Ejecutar CV completo
  result <- run_cox_cv(data, formula = vars)
  
  # Verificaciones básicas de estructura
  expected_names <- c("fold_results", "global_c_index", "mean_c_index", 
                     "sd_c_index", "n_folds", "formula")
  stopifnot(all(expected_names %in% names(result)))
  
  # Verificar fold_results
  stopifnot(is.list(result$fold_results))
  stopifnot(length(result$fold_results) == 4)  # 4 folds
  
  # Verificar métricas globales
  stopifnot(is.numeric(result$global_c_index))
  stopifnot(is.numeric(result$mean_c_index))
  stopifnot(is.numeric(result$sd_c_index))
  stopifnot(result$n_folds == 4)
  
  # Verificar rangos
  stopifnot(result$global_c_index >= 0 && result$global_c_index <= 1)
  stopifnot(result$mean_c_index >= 0 && result$mean_c_index <= 1)
  stopifnot(result$sd_c_index >= 0)
  
  # Verificar que cada fold tiene resultados
  for (fold_result in result$fold_results) {
    stopifnot("c_index" %in% names(fold_result))
    stopifnot(is.numeric(fold_result$c_index))
  }
  
  cat("  ✓ CV completo ejecuta sin errores\n")
  cat("  ✓ Estructura de output correcta\n")
  cat("  ✓ Métricas en rangos válidos\n")
  cat("  ✓ Resultados por fold incluidos\n")
}

# Test 4: Verificar manejo de diferentes formatos de formula
test_formula_formats <- function() {
  cat("\nTest 4: Formatos de formula - vector vs formula\n")
  
  data <- create_synthetic_data()
  
  # Test con vector de variables
  vars_vector <- c("AGE", "HB", "TP53")
  result1 <- train_cox_fold(data, fold_id = 1, formula = vars_vector)
  
  # Test con formula explícita
  formula_obj <- as.formula("Surv(LFS_YEARS, LFS_STATUS) ~ AGE + HB + TP53")
  result2 <- train_cox_fold(data, fold_id = 1, formula = formula_obj)
  
  # Ambos deben funcionar
  stopifnot(inherits(result1$model, "coxph"))
  stopifnot(inherits(result2$model, "coxph"))
  
  # Los coeficientes deben ser similares (mismo modelo)
  coef1 <- coef(result1$model)
  coef2 <- coef(result2$model)
  stopifnot(length(coef1) == length(coef2))
  
  cat("  ✓ Vector de variables funciona\n")
  cat("  ✓ Formula explícita funciona\n")
  cat("  ✓ Ambos formatos producen modelos válidos\n")
}

# Test 5: Verificar manejo de casos edge
test_edge_cases <- function() {
  cat("\nTest 5: Casos extremos\n")
  
  # Datos muy pequeños
  small_data <- create_synthetic_data(n = 15, n_folds = 3)
  
  # Debe funcionar incluso con pocos datos
  result <- run_cox_cv(small_data, formula = c("AGE", "TP53"))
  
  stopifnot(is.list(result))
  stopifnot("global_c_index" %in% names(result))
  
  # Verificar que no hay NAs en métricas principales
  stopifnot(!is.na(result$global_c_index))
  stopifnot(!is.na(result$mean_c_index))
  
  cat("  ✓ Funciona con datos pequeños\n")
  cat("  ✓ No produce NAs en métricas principales\n")
}

# Función principal para ejecutar todos los tests
run_all_tests <- function() {
  cat("=== EJECUTANDO TESTS UNITARIOS - COX CV FUNCTIONS ===\n\n")
  
  # Verificar que survival está disponible
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' no está disponible")
  }
  
  tryCatch({
    test_train_cox_fold_basic()
    test_evaluate_cox_fold_basic()
    test_run_cox_cv_basic()
    test_formula_formats()
    test_edge_cases()
    
    cat("\n=== TODOS LOS TESTS PASARON EXITOSAMENTE ===\n")
    cat("✅ Las funciones Cox CV están funcionando correctamente\n")
    
  }, error = function(e) {
    cat("\n❌ ERROR EN TESTS:\n")
    cat(paste("Error:", e$message, "\n"))
    stop("Tests fallaron")
  })
}

# Ejecutar tests si se corre directamente
if (!interactive()) {
  run_all_tests()
}
