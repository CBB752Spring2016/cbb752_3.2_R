# CBB752 Final Project 3.2, R card, by Julian Q Zhou

## Objective

Given a tab-delimited protein-protein interaction file of MITAB2.5 format, calculate degree centrality and betweenness centrality for each unique interactor/protein, and export the output as a csv file.

## Source code

Available [here](https://github.com/jqz752/cbb752_3.2_R)

* `r_centrality.R`: main script

## Sample input
* A MITAB2.5 file that looks like below (e.g. `sample_input_Rnorv20160114.txt`). It is assumed that both `ID interactor A` and `ID interactor B` contain DIP identifiers of the format `DIP-xN`, where `x` is a number.

|ID interactor A	| ID interactor B	| Alt. ID interactor A	| Alt. ID interactor B	| Alias(es) interactor A	| Alias(es) interactor B	| Interaction detection method(s)	| Publication 1st author(s)	| Publication Identifier(s)	| Taxid interactor A	| Taxid interactor B	| Interaction type(s)	| Source database(s)	| Interaction identifier(s)	| Confidence value(s)	| Processing Status	|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|DIP-383N\|refseq:NP_001094\|uniprotkb:P35609	| DIP-674N\|refseq:NP_058706\|uniprotkb:P35439	| -|	-|	-|	-|	MI:0019(coimmunoprecipitation)|	-	|pubmed:9009191\|pubmed:DIP-290S	|taxid:9606(Homo sapiens)	|taxid:10116(Rattus norvegicus)	|MI:0218(physical interaction)|	MI:0465(dip)|	DIP-30E	|dip-quality-status:core|	dip:0002(small scale)	|	-|
|DIP-1072N\|refseq:NP_059040\|uniprotkb:P18265	| DIP-736N\|refseq:NP_033782\|uniprotkb:P31750	| -	|-	|-	|-	|MI:0045(experimental interaction detection)|	-	|pubmed:9005851\|pubmed:DIP-168S|	taxid:10116(Rattus norvegicus)|	taxid:10090(Mus musculus)	|MI:0218(physical interaction)|	MI:0465(dip)|	DIP-93E|	dip-quality-status:core	|dip:0002(small scale)		-

## Sample output
* A comma-separated .csv file that looks like below (e.g. `sample_output.csv`):

interactor|degree_centrality|degree_centrality_norm|betweenness_centrality|betweenness_centrality_norm
---|---|---|---|---
DIP-383N|1|0.0015600624024961|0|0
DIP-358N|4|0.0062402496099844|32.6666666666667|0.00015925637025481
DIP-5717N|8|0.0124804992199688|643|0.0031347503900156
DIP-5790N|2|0.0031201248049922|1|4.87519500780031e-06

## Usage

`get.centrality(input.ppi, output.csv)`

* `input.ppi`: filename of input ppi file; must be a tab-delimited file in MITAB2.5 format.
* `output.csv`: filename of output csv file

Only non-self interactions is considered when constructing the adjacency matrix. That is, an interaction in which `ID interactor A` and `ID interactor B` are the same does not form an edge, and therefore does not contribute to that interactor/protein's degree centrality or betweenness centrality.

Degree centrality of an interactor is normalized by `n-1`, and betweenness centrality of an interactor is normalized by `(n-1)(n-2)/2`, where `n` is the number of unique interactors.

Betweenness centrality is calculated using Brandes algorithm ([Brandes, 2001](http://algo.uni-konstanz.de/publications/b-fabc-01.pdf)), implemented via `brandes.betweenness.centrality()` from the `RBGL` package ([Carey, Long & Gentleman, 2016](https://www.bioconductor.org/packages/release/bioc/html/RBGL.html)). When running `get.centrality`, it checks for presence of `RBGL` and installs it if absent before calling `brandes.betweenness.centrality()`. `brandes.betweenness.centrality()` takes in an object of the class `graph`, converted from the adjacency matrix constructued, and computes absolute betweenness centrality using Brandes algorithm.

#### Example
`> get.centrality(input.ppi = "sample_input_Rnorv20160114.txt", output.csv = "sample_output.csv")`

`[1] "reading in sample_input_Rnorv20160114.txt..."`

`[1] "constructing adjacency matrix for 642 unique interactors (this might take a while)..."`

`[1] "computing degree centrality for each unique interactor..."`

`[1] "computing betweenness centrality for each unique interactor..."`

`Loading required package: RBGL`

`Loading required package: graph`

`[1] "exporting..."`
