# Using the normalized, non-log transformed data
df_expr_matrix = integrated.X # integrated.X ==> 지금 sparse matrix가 아니라 numpy.matrix로 되어있음.
df_expr_matrix = df_expr_matrix.T
df_expr_matrix = pd.DataFrame(df_expr_matrix)
df_expr_matrix.columns = test3.obs.index
df_expr_matrix.set_index(test3.var.index, inplace=True)

