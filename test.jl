module test
x = "abc def"

x = "1.234"
x = parse(Float64, x)
N = Int(1e2)

a = [[1,2,3],[4,5,6],[7,8,9]]
for ai = a
    @show ai
end

end