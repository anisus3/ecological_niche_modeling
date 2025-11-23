## ==========================================================
##  AUC (ROC) для BIOMOD моделей + сохранение в CSV
##  объект с моделями: models  (BIOMOD.models.out)
## ==========================================================

library(biomod2)
library(dplyr)

# --- 0. Папка для вывода ---
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1. Извлекаем оценки моделей ---
# В твоей версии biomod2 get_evaluations() -> data.frame
evals <- get_evaluations(models)

# При желании можно посмотреть структуру:
# str(evals)
# head(evals)

# Ожидаемые ключевые столбцы:
# full.name, PA, run, algo, metric.eval,
# cutoff, sensitivity, specificity,
# calibration, validation, evaluation

# --- 2. Фильтруем только AUC (ROC) ---
# В biomod2 AUC обычно обозначается "ROC"
auc_runs <- evals %>%
  dplyr::filter(metric.eval %in% c("ROC", "AUC", "AUCroc"))

# --- 3. Полная таблица для supplementary (все прогоны) ---
auc_all_out <- auc_runs %>%
  dplyr::select(
    full.name,      # уникальное имя модель+прогон
    PA,             # набор псевдоотрицаний (если есть)
    run,            # номер прогона
    algo,           # алгоритм (GLM, RF, GBM, MAXENT...)
    metric.eval,    # ROC / AUC
    cutoff, sensitivity, specificity,
    calibration,    # AUC на обучающих
    validation,     # AUC на CV
    evaluation      # AUC на независимых данных (если задавал)
  )

write.csv(
  auc_all_out,
  file = file.path(out_dir, "AUC_all_runs_models.csv"),
  row.names = FALSE
)

cat("Сохранено: ", file.path(out_dir, "AUC_all_runs_models.csv"), "\n")

# --- 4. Сводка по алгоритмам (для текста статьи) ---
# Чаще всего в тексте дают средний AUC по cross-validation (validation)

auc_summary_algo <- auc_runs %>%
  group_by(algo) %>%
  summarise(
    N_runs          = sum(!is.na(validation)),
    AUC_calib_mean  = mean(calibration, na.rm = TRUE),
    AUC_calib_sd    = sd(calibration, na.rm = TRUE),
    AUC_valid_mean  = mean(validation, na.rm = TRUE),
    AUC_valid_sd    = sd(validation, na.rm = TRUE),
    AUC_eval_mean   = mean(evaluation, na.rm = TRUE),
    AUC_eval_sd     = sd(evaluation, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  arrange(desc(AUC_valid_mean))

write.csv(
  auc_summary_algo,
  file = file.path(out_dir, "AUC_summary_by_algo.csv"),
  row.names = FALSE
)

cat("Сохранено: ", file.path(out_dir, "AUC_summary_by_algo.csv"), "\n")

# --- 5. (Опционально) Сводка по отдельным моделям (full.name) ---
auc_summary_model <- auc_runs %>%
  dplyr::select(
    full.name, PA, run, algo,
    AUC_calibration = calibration,
    AUC_validation  = validation,
    AUC_evaluation  = evaluation
  )

write.csv(
  auc_summary_model,
  file = file.path(out_dir, "AUC_summary_by_model.csv"),
  row.names = FALSE
)

cat("Сохранено: ", file.path(out_dir, "AUC_summary_by_model.csv"), "\n")
