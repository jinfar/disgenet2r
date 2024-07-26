# disgenet2r

## This is fork of an old version, new versions can be found here 
## https://gitlab.com/medbio/disgenet2r

`disgenet2r` is an R package to query and expand DisGeNET data (www.disgenet.org), and to visualize the results within R framework.
The disgenet2r is designed to query data for DisGeNET v7.0 (May, 2020).
## What is this repository for?

This report is used for package distribution and testing until it is ready to be published in BioConductor.

## Package' Status

 * __Version__: 0.99.3
 * __Subversion__: `20230223`
 * __Authors__:  IBI group
 * __Maintainer__: <support@disgenet.org>

## How to start

### Installation

The package, `disgenet2r` can be installed using `devtools` from this repository:

```R
library(devtools)
install_bitbucket("ibi_group/disgenet2r")
```

### Obtaining the API key

Before using the package, you will need to create a free DisGeNET account at [http://disgenet.org/signup](http://disgenet.org/signup). Once you have completed the registration process, use the **get_disgenet_api_key** function to retrieve your API key.


```R
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
                  email = "user@gmail.com", 
                  password = "myspwd" )
```
After retrieving the API key, run the line below so the key is available for all the **disgenet2r** functions. 

```R
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)
```


### Querying DisGeNET:

The following lines show two examples of how DisGeNET can be queried using `disgenet2r`:

 * __Gene Query__

```R
library(disgenet2r)
gq <- gene2disease(gene = 3953, 
 vocabulary = "ENTREZ",
    database = "ALL", 
    score = c( 0.1,1)
)
```

 * __Disease Query__

```R
library(disgenet2r)
dq <- disease2gene(disease = "C0028754", 
    database = "ALL",
    score = c(0.3,1) 
)
```

A detailed documentation of the functions of the package is available at http://www.disgenet.org/disgenet2r.

## COPYRIGHT

Copyright (C) 2023 IBI group.

disgenet2r is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

disgenet2r is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

## How to cite disgenet2r:

Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno, E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for disease genomics: 2019 update. Nucleic acids research, 48(D1), D845-D855.


