load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")



cc_binary(
    name = "exe",
    srcs = glob([
        "src/*.cpp",
        "src/*.h",
        "src/*.hpp",
    ], allow_empty = True),
    deps = [
        "@sfml_src//:sfml",
        "@glm_src//:glm",
        "@imgui_src//:imgui",
        "@imgui_sfml_src//:imgui_sfml",
        "@tracy_src//:tracy",
        "@nlohmann_json//:json",
    ],
    defines = ["SFML_STATIC"],
    copts = [
        "-std=c++26",
        "-fopenmp",
    ],
    linkopts = [
        "-fopenmp",
        "-lX11",
        "-lXrandr", 
        "-lXi",
        "-lXcursor",
        "-lxcb",
        "-lxcb-image",
        "-ludev",
        "-lGL",
        "-lfreetype",
        "-lpthread",
        "-ldl",
    ],
)