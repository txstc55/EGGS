python do_test.py --e 0 --d 1 --o slim_demo1 
python do_test.py --e 0 --d 1 --o slim_demo1_400k_triangles --f face_400k.obj 
python do_test.py --e 0 --d 2 --o slim_demo2 
python do_test.py --e 0 --d 3 --o slim_demo3 
python do_test.py --e 1 --w 75 --o cot_smoothing_10k_surface_37385 --f 37385_sf.obj 
python do_test.py --e 1 --o cot_smoothing_armadillo 
python do_test.py --e 2 --d 0 --o optical_flow_small_example 
python do_test.py --e 2 --d 1 --o optical_flow_large_example
python do_test.py --e 3 --o cot_matrix_armadillo 
python do_test.py --e 3 --o cot_matrix_37385_sf --f 37385_sf.obj 

python do_numeric_test.py
