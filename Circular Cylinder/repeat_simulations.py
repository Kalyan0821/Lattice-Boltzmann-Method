import os

for N0 in [50, 100, 200, 300, 400]:
    u0 = 0.15
    nul = 3*N0/400  # maintain Re=20
    assert u0*N0/nul == 20

    for obstacle_mode in ["bb", "ns", "curved"]:
        base = "./output-Re{20}_u{" + str(u0) + "}_N0{" + str(N0) + "}_nu{" + str(nul) + "}_"

        if obstacle_mode == "curved":
            for curved_scheme in [1, 2]:
                path = base + "C" + str(curved_scheme)
                if not os.path.exists(path):
                    os.makedirs(path)
                    os.system("g++ -std=c++17 -fopenmp CircularCylinder.cpp StairBoundary.cpp CurvedBoundary.cpp -o out")
                    os.system(f"./out {u0} {N0} {nul} {obstacle_mode} {curved_scheme}")
        else:
            path = base + obstacle_mode
            if not os.path.exists(path):
                os.makedirs(path)
                os.system("g++ -std=c++17 -fopenmp CircularCylinder.cpp StairBoundary.cpp CurvedBoundary.cpp -o out")
                os.system(f"./out {u0} {N0} {nul} {obstacle_mode} {0}")

