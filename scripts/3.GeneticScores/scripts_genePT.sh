## Full input dataset
### Include all genes, Cyto IPSSR, individual karyotypes and all clinical variables (exclude MONOCYTES and PB_BLASTS due to missings)
docker run --rm --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.7 bash

python scripts/3.GeneticScores/prepare_data_unified.py \
    --config scripts/3.GeneticScores/configs/full_model_genePT.yaml \
    --output results/gnn/preprocess/full_input_genePT/


python scripts/3.GeneticScores/run_gnn_cv.py \
    --graph_data_path results/gnn/preprocess/full_input_genePT/full_model.pt \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --results_dir results/gnn/full_input_model_genePT \
    --lightning_logs_dir lightning_logs/full_input_model_genePT


python scripts/3.GeneticScores/evaluate_cv_results.py --results_dir results/gnn/full_input_model_genePT \
    --output_dir results/gnn/evaluations/full_model_genePT \
    --external_benchmarks IPSSM:results/gnn/IPSSM_cv_folds/cv_summary.tsv \
    --verbose

python scripts/3.GeneticScores/apply_GNN_models_cv.py \
    --graph_data_path results/gnn/preprocess/full_input_genePT/full_model.pt \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --check_point_path results/gnn/full_input_model_genePT/model_checkpoints/ \
    --results_dir results/gnn/full_input_model_genePT/ \
    --models v1_boolean v1_vaf v2_boolean v2_vaf v3_boolean v3_vaf v4_boolean v4_vaf

# Create graphs per individual
python scripts/3.GeneticScores/define_mutation_combinations.py \
    --config scripts/3.GeneticScores/configs/full_model_inference_genePT.yaml \
    --output results/gnn/preprocess/full_input_genePT/

## Define reference scores for each model
python scripts/3.GeneticScores/define_reference_score.py \
    --graph_data_path results/gnn/preprocess/full_input_genePT/full_model.pt \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --check_point_path results/gnn/full_input_model_genePT/model_checkpoints/ \
    --model_name v1_boolean \
    --fold 0 

python scripts/3.GeneticScores/evaluate_gene_effect.py \
    --graph_data_path results/gnn/preprocess/full_input_genePT/full_model_individual.pt \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --check_point_path results/gnn/full_input_model_genePT/model_checkpoints/fold_0/v1_boolean_fold_0.ckpt \
    --model_name v1_boolean \
    --reference_score -1.14249849319458 \
    --results_dir results/gnn/full_input_model_genePT/v1_boolean/

## Define graphs by gene
python scripts/3.GeneticScores/define_graphs_by_gene.py \
    --config scripts/3.GeneticScores/configs/full_model_inference_genePT.yaml \
    --output results/gnn/preprocess/genePT_per_gene/


python scripts/3.GeneticScores/run_gnn_genes.py \
    --graph_data_path results/gnn/preprocess/genePT_per_gene/full_model_individual.pt \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --results_dir results/gnn/genePT_model_per_gene \
    --models v1_boolean \
    --lightning_logs_dir lightning_logs/genePT_model_per_gene/

python scripts/3.GeneticScores/evaluate_gene_effect_per_gene.py \
    --graph_data_path results/gnn/preprocess/full_input_genePT/full_model_individual.pt \
    --graph_gene_path results/gnn/preprocess/genePT_per_gene/full_model_individual.pt  \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --check_point_path results/gnn/genePT_model_per_gene/model_checkpoints/ \
    --model_name v1_boolean \
    --results_dir results/gnn/genePT_model_per_gene/v1_boolean/