#!/bin/bash
# Script de ejemplo para ejecutar el pipeline con Docker
# 
# Este script demuestra diferentes formas de usar el pipeline MDS con Docker

set -e  # Salir si cualquier comando falla

# Variables de configuración
DOCKER_IMAGE="mds_subtypes_python:1.6"
WORKSPACE_DIR="/home/cruiza/data/cruiza/MDS_subtypes"
SCRIPTS_DIR="/workspace/scripts/3.GeneticScores"

echo "🐳 Pipeline MDS con Docker - Ejemplos de Uso"
echo "=============================================="

# Función para verificar Docker
check_docker() {
    if ! command -v docker &> /dev/null; then
        echo "❌ Error: Docker no está instalado o no está en el PATH"
        exit 1
    fi
    
    if ! docker images | grep -q "$DOCKER_IMAGE"; then
        echo "❌ Error: Imagen Docker no encontrada: $DOCKER_IMAGE"
        echo "💡 Construye la imagen con: docker build -t $DOCKER_IMAGE ."
        exit 1
    fi
    
    echo "✅ Docker y la imagen están disponibles"
}

# Función para ejecutar el pipeline básico
run_basic_pipeline() {
    echo ""
    echo "📊 Ejecutando pipeline básico..."
    echo "--------------------------------"
    
    docker run --rm \
        -v "$WORKSPACE_DIR:/workspace" \
        "$DOCKER_IMAGE" \
        python "$SCRIPTS_DIR/prepare_data_unified.py" \
        --config "$SCRIPTS_DIR/configs/real_data_config.yaml" \
        --output "/workspace/results/docker_basic" \
        --verbose
    
    echo "✅ Pipeline básico completado"
}

# Función para ejecutar validación de calidad
run_quality_validation() {
    echo ""
    echo "🔍 Ejecutando validación de calidad..."
    echo "------------------------------------"
    
    docker run --rm \
        -v "$WORKSPACE_DIR:/workspace" \
        "$DOCKER_IMAGE" \
        python "$SCRIPTS_DIR/tests/validate_step3c.py" \
        --data-dir "/workspace/results/docker_basic" \
        --verbose
    
    echo "✅ Validación de calidad completada"
}

# Función para ejecutar validación de compatibilidad
run_compatibility_validation() {
    echo ""
    echo "🧪 Ejecutando validación de compatibilidad..."
    echo "--------------------------------------------"
    
    docker run --rm \
        -v "$WORKSPACE_DIR:/workspace" \
        "$DOCKER_IMAGE" \
        python "$SCRIPTS_DIR/tests/validate_step3d.py" \
        --data-dir "/workspace/results/docker_basic" \
        --verbose
    
    echo "✅ Validación de compatibilidad completada"
}

# Función para ejecutar análisis completo automatizado
run_full_automated_analysis() {
    echo ""
    echo "🚀 Ejecutando análisis completo automatizado..."
    echo "----------------------------------------------"
    
    docker run --rm \
        -v "$WORKSPACE_DIR:/workspace" \
        "$DOCKER_IMAGE" \
        python "$SCRIPTS_DIR/tests/run_pipeline_and_validate.py" \
        --config "$SCRIPTS_DIR/configs/real_data_config.yaml" \
        --output "/workspace/results/docker_automated" \
        --run-quality-validation \
        --run-model-validation
    
    echo "✅ Análisis completo automatizado completado"
}

# Función para mostrar contenido de un contenedor interactivo
interactive_session() {
    echo ""
    echo "💻 Iniciando sesión interactiva..."
    echo "--------------------------------"
    echo "Usa 'exit' para salir del contenedor"
    
    docker run -it --rm \
        -v "$WORKSPACE_DIR:/workspace" \
        "$DOCKER_IMAGE" \
        bash
}

# Función para mostrar ayuda
show_help() {
    echo ""
    echo "📖 Uso: $0 [comando]"
    echo ""
    echo "Comandos disponibles:"
    echo "  basic       - Ejecutar pipeline básico"
    echo "  quality     - Validar calidad de datos (requiere datos previos)"
    echo "  compat      - Validar compatibilidad con modelos (requiere datos previos)"
    echo "  full        - Ejecutar análisis completo automatizado"
    echo "  interactive - Iniciar sesión interactiva en contenedor"
    echo "  help        - Mostrar esta ayuda"
    echo ""
    echo "Ejemplos:"
    echo "  $0 basic"
    echo "  $0 full"
    echo "  $0 interactive"
}

# Función principal
main() {
    local command="${1:-help}"
    
    # Verificar Docker
    check_docker
    
    case "$command" in
        "basic")
            run_basic_pipeline
            echo ""
            echo "💡 Para validar los resultados, ejecuta:"
            echo "   $0 quality"
            echo "   $0 compat"
            ;;
        "quality")
            if [ ! -d "$WORKSPACE_DIR/results/docker_basic" ]; then
                echo "❌ Error: No se encontraron datos para validar"
                echo "💡 Ejecuta primero: $0 basic"
                exit 1
            fi
            run_quality_validation
            ;;
        "compat")
            if [ ! -d "$WORKSPACE_DIR/results/docker_basic" ]; then
                echo "❌ Error: No se encontraron datos para validar"
                echo "💡 Ejecuta primero: $0 basic"
                exit 1
            fi
            run_compatibility_validation
            ;;
        "full")
            run_full_automated_analysis
            echo ""
            echo "🎉 Análisis completo finalizado!"
            echo "📁 Revisa los resultados en: $WORKSPACE_DIR/results/docker_automated/"
            ;;
        "interactive")
            interactive_session
            ;;
        "help"|"--help"|"-h")
            show_help
            ;;
        *)
            echo "❌ Comando desconocido: $command"
            show_help
            exit 1
            ;;
    esac
}

# Ejecutar función principal con todos los argumentos
main "$@"
