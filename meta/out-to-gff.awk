BEGIN { print "##gff-version 3"; }
{
    s = substr($1, 1, 1)
    if (s == ">") {
        seqid = substr($1, 2)
    } else {
        s = split($0, t, "\t")
        id = "ID=" seqid "_" t[1] "_" t[2] "_" t[3] ";product=predicted protein"
        print seqid "\tFGS\tCDS\t" t[1] "\t" t[2] "\t.\t" t[3] "\t" int(t[4] - 1) "\t" id
    }
}