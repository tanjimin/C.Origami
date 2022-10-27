import sys
import torch
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import corigami.model.corigami_models as corigami_models

def load_default(model_path):
    model_name = 'ConvTransModel'
    mid_hidden = 256
    model = get_model(model_name, mid_hidden)
    load_checkpoint(model, model_path)
    return model

def get_model(model_name, mid_hidden, num_genomic_features=2):
    ModelClass = getattr(corigami_models, model_name)
    model = ModelClass(num_genomic_features, mid_hidden = mid_hidden)
    return model

def load_checkpoint(model, model_path):
    print('Loading weights')
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    checkpoint = torch.load(model_path, map_location=device)
    model_weights = checkpoint['state_dict']

    # Edit keys
    for key in list(model_weights):
        model_weights[key.replace('model.', '')] = model_weights.pop(key)
    model.load_state_dict(model_weights)
    model.eval()
    return model

if __name__ == '__main__':
    main()
