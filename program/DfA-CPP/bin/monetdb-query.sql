select * from data_transformation;

SELECT * FROM data_set;

SELECT * FROM data_dependency;

SELECT * FROM extractor;

SELECT * FROM attribute;

SELECT a.name, a.type, s.tag, a.extractor_id FROM attribute a, data_set s WHERE a.ds_id=s.id ;