%% script to test_PointProcessData

D = LoadDataIntoProcess(1);
maxlags = 500;
D.Data.PreProcessCP( maxlags );
