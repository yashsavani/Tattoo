# Tattoo Generator

This is a personal tattoo generator project. We use a Delaunay triangulation of five randomly 
sampled points from a Dirichlet distribution. We also manage the randomness by selecting the 
top distances that are calculated from the volumes of the triangles and a huber loss on
the angles of the triangles.

### Prerequisites

All the libraries needed to run the project are provided in the requirements.txt file

```
pip install requirements.txt
```

## Getting Started

To run the code simply run
```
python tattoo.py
```
The results will be saved in the `figs` directory.

## Authors

* **Yash Savani** - *Initial work* - [Yash Savani](https://github.com/yashsavani)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
