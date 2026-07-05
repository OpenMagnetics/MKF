import pandas
import ndjson
import os

from sqlalchemy import create_engine
from sqlalchemy.engine import URL
from sqlalchemy import Column, Integer, String, DateTime, JSON, Float
from sqlalchemy.orm import declarative_base
from datetime import datetime
import sqlalchemy
import pandas
from sqlalchemy.orm import sessionmaker

current_advanced_core_materials_path = "/home/alf/OpenMagnetics/MAS/data/advanced_core_materials.ndjson"

current_advanced_core_materials = pandas.DataFrame(ndjson.load(open(current_advanced_core_materials_path, "r")))
current_advanced_core_materials = current_advanced_core_materials.where(pandas.notnull(current_advanced_core_materials), None)
current_advanced_core_materials["createdAt"] = [datetime.now()] * len(current_advanced_core_materials.index)
current_advanced_core_materials["updatedAt"] = [datetime.now()] * len(current_advanced_core_materials.index)

Base = declarative_base()

driver = "postgresql"
address = os.getenv('OM_DB_ADDRESS')
port = os.getenv('OM_DB_PORT')
name = os.getenv('OM_DB_NAME')
user = os.getenv('OM_DB_USER')
password = os.getenv('OM_DB_PASSWORD')

engine = sqlalchemy.create_engine(f"{driver}://{user}:{password}@{address}:{port}/{name}")

current_advanced_core_materials.to_sql(
    'advanced_core_materials',
    con=engine,
    if_exists='replace',
    index=True,
    dtype={
        "manufacturerInfo": sqlalchemy.types.JSON, 
        "permeability": sqlalchemy.types.JSON, 
        "volumetricLosses": sqlalchemy.types.JSON, 
        "bhCycle": sqlalchemy.types.JSON, 
        "massLosses": sqlalchemy.types.JSON, 
    }
)  

engine.execute('ALTER TABLE public.advanced_core_materials ADD PRIMARY KEY (index);')
