mainList = Channel.from('314','246','234')
    .toList()


refList = mainList.toType(String)//.join('","')
//refList = '[' + mainList.join(',') + ']'


process sayHello {
executor = 'local'
    """
    printf 'Hello world! ${refList} \n'
    """
}