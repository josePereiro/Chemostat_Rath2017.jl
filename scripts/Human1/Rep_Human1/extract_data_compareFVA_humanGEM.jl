# ### comparativeFVA
# ---

# +
# This is a function from GECKO package (https://github.com/SysBioChalmers/GECKO.git).
# comparativeFVA(model,ecModel,CsourceUptk,false,1E-8);
# where model is the tINIT_orig model, ecModel is the model with enzymatic constraints, 
# CsourceUptk = 'HMR_9034', the 'false' parameter refers to 'chemostat' condition and 1-E8 to 
# 'tol' (numerical tolerance for a flux and variability range to be considered as zero)
# -

# ### Orig model

# # Gets the optimal value for ecirrevModel and fixes the objective value to
# # this for both models
# println("\nFixing objective ($obj_ider)")
# obj_val = Ch.LP.fba(ec_model, obj_ider).obj_val
# println("\tobj_val: ", obj_val)
# for var in [:orig_model, :ec_model]
#     model = eval(var)comp
#     Ch.Utils.bounds!(model, obj_ider, obj_val - zeroth, obj_val + zeroth)
# end

# println("\n Check if models are feasible")
# for var in [:orig_model, :ec_model]
#     println(var)
#     model = eval(var)
#     obj_val = Ch.LP.fba(model, obj_ider).obj_val;
#     if obj_val < zeroth
#         println("\tWARNING: Constrained $var is unfeasible")
#     else
#         println("\tConstrained $var is feasible, objval: ", obj_val)
#     end
# end



# +
# # TODO: package this
# fl = "orig_model_fva.bson"
# tagsave(fl, Dict("dat" => (lvals, uvals)))
# println(relpath(fl), " created, size: ", filesize(fl), " bytes")

# +
# diff_ = (uvals .- lvals);
# nz_diff_ = diff_[abs.(diff_) .> zeroth];
# cdist = ecdf(nz_diff_);
# xs = 10.0.^collect(-8:0.5:3.5)
# ys = cdist.(xs)

# p = Plots.plot(xscale = :log10)
# p = Plots.plot!(p, xs, ys, label = "")
# p = Plots.scatter!(p, [mean(diff_)],[0.5], label = "", ms = 10)
# -

# mean(diff_)


