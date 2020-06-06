python plot.py result_datas/slim_demo1.json
python plot.py result_datas/slim_demo1_400k_triangles.json
python plot.py result_datas/slim_demo2.json
python plot.py result_datas/slim_demo3.json
python plot.py result_datas/cot_smoothing_armadillo.json
python plot.py result_datas/cot_smoothing_10k_surface_37385.json
python plot.py result_datas/optical_flow_small_example.json
python plot.py result_datas/optical_flow_large_example.json
python plot.py result_datas/cot_matrix_armadillo.json
python plot.py result_datas/cot_matrix_37385_sf.json

python plot_numeric.py

python plot_opt.py result_datas/optical_flow_small_example_u.txt result_datas/optical_flow_small_example_v.txt tutorial/data/box.0.bmp tutorial/data/box.1.bmp test_result_graphs/implemented_small_opt.pdf
python plot_opt.py result_datas/optical_flow_large_example_u.txt result_datas/optical_flow_large_example_v.txt tutorial/data/frame10.png tutorial/data/frame11.png test_result_graphs/implemented_large_opt.pdf