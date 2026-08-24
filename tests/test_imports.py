import importlib


def test_process_experiment_module_imports():
    module = importlib.import_module('process_experiment')
    assert hasattr(module, 'process_experiment')
    assert hasattr(module, 'main')
