# I need the following so that during "check" R doesnt complain about index in lgc
# and for Random.Seed in set_rand_state
utils::globalVariables(".Random.seed")
utils::globalVariables("index")
