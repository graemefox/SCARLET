from NN_model import NN_classifier
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Args for classify_NN_bedmethyl')
parser.add_argument('-m','--model', \
                    help='The model to use', \
                    required=True)
parser.add_argument('-i','--input_bed', \
                    help='The input methylBed', \
                    required=True)
parser.add_argument('-v', '--votes', \
                    help='The votes file to write', \
                    required=True)
parser.add_argument('-o', '--output', \
                    help='Output file', \
                    required=True)
args = parser.parse_args()

# give args snappy names
model = args.model
bed = args.input_bed
votes = args.votes
output = args.output

#NN = NN_classifier(snakemake.input["model"])
NN = NN_classifier(model)

#predictions, class_labels, n_features = NN.predict_from_bedMethyl(snakemake.input["bed"])
predictions, class_labels, n_features = NN.predict_from_bedMethyl(bed)

# write predictions to table
df = pd.DataFrame({'class': class_labels, 'score': predictions, 'num_features': n_features})
#df.to_csv(snakemake.output['votes'], sep='\t')
df.to_csv(votes, sep='\t')

print(df)

# write summary to txt file
summary = ['Number of features: ' + str(n_features),
           'Predicted Class: ' + str(class_labels[0]),
           "Score: " + str(predictions[0])
          ]

#with open(snakemake.output["txt"], 'w') as f:
with open(output, 'w') as f:
  f.write("\n".join(summary))
