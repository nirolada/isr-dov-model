from typing import get_type_hints


def validate_input(obj, **kwargs):
    hints = get_type_hints(obj)
    for attr_name, attr_type in hints.items():
        if attr_name == 'return':
            continue
        if not isinstance(kwargs[attr_name], attr_type):
            raise TypeError(f'Argument {attr_name} must be {attr_type}')